# Gera projeto_final_semsaida.ipynb (sem saídas). Execute depois com run_nb.py
import os, sys, nbformat
from nbformat.v4 import new_notebook, new_markdown_cell as md, new_code_cell as code

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
PLANO_B = "(será preenchido após o primeiro download)"
fasta = os.path.join(RAIZ, "dados", "cytb_especies.fasta")
if os.path.exists(fasta):
    from Bio import SeqIO
    escolhidas = ["Humano", "Chimpanzé", "Gorila", "Cachorro", "Gato", "Vaca", "Baleia-azul", "Galinha", "Peixe-zebra"]
    linhas = []
    for r in SeqIO.parse(fasta, "fasta"):
        nome = r.description.split(" | ")[0].split(" ", 1)[1]
        if nome in escolhidas:
            linhas.append(f">{r.description}")
            s = str(r.seq)
            linhas += [s[i:i+80] for i in range(0, len(s), 80)]
    PLANO_B = "\n".join(linhas)

C = []

C.append(code("""# Instalação do Biopython (no Colab ele não vem instalado). Rode esta célula primeiro.
!pip install -q biopython"""))

C.append(md("""# Quem é parente de quem? A árvore da vida pelo DNA

**Ambientes em Computação — UFMG · Projeto final**

Neste notebook comparamos o gene mitocondrial **citocromo b (CYTB)** de 15 espécies, usando um algoritmo de alinhamento
**implementado do zero** (Needleman-Wunsch), e a partir dele construímos uma **matriz de distância**, um **heatmap** e uma **árvore filogenética**.

**Roteiro do notebook**
1. Introdução
2. Dados (download do NCBI)
3. Algoritmo (Needleman-Wunsch)
4. Análise (distâncias, heatmap, árvore)
5. Complexidade (tempo × tamanho)
6. Conclusões e limitações"""))

# ---------------- 1. Introdução
C.append(md("""## 1. Introdução

### A pergunta central
> **Dá para descobrir quem é parente de quem apenas comparando o DNA das espécies?**

Se dois animais têm um ancestral comum recente, seus DNAs tiveram pouco tempo para acumular diferenças. Quanto mais distante o parentesco, mais diferenças acumuladas.
Então, *medir o quanto dois DNAs se parecem* é uma forma de medir parentesco.

### DNA como um "texto"
O DNA pode ser visto como um texto escrito com apenas **4 letras**: **A, C, G, T**. Um gene é como uma frase longa nesse alfabeto.
O gene que usaremos, o **CYTB** (citocromo b), fica nas mitocôndrias (as "usinas de energia" das células) e tem cerca de **1.140 letras**.
Ele é muito usado em estudos de parentesco porque todos os vertebrados o possuem, e ele muda devagar o bastante para ser comparável.

### O que é "alinhar" duas sequências?
Imagine duas frases quase iguais:

```
O GATO COMEU O PEIXE
O GATO COMEU - PEIXE      <- faltou uma palavra; marcamos o "buraco" com '-'
```

**Alinhar** é colocar as duas sequências lado a lado, inserindo "buracos" (*gaps*) quando necessário, de modo que as partes parecidas fiquem na mesma coluna.
Em seguida, contamos quantas colunas são iguais (**identidade**). No DNA, um buraco representa uma inserção ou remoção de letras ao longo da evolução, e uma coluna diferente representa uma mutação.

O algoritmo de **Needleman-Wunsch** encontra o *melhor* alinhamento possível de ponta a ponta (alinhamento **global**), usando uma pontuação:
**+ para letras iguais (match), − para letras diferentes (mismatch) e − para buracos (gap)**."""))

# ---------------- Setup
C.append(code("""# ===== Importações e configurações gerais =====
import os, io, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from Bio import Entrez, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform
from IPython.display import display

# >>> EDITE AQUI: o NCBI pede um e-mail de contato para quem usa a API <<<
ENTREZ_EMAIL = "seu_email@exemplo.com"
Entrez.email = ENTREZ_EMAIL

# Pastas de saída (no Colab aparecem no painel de arquivos à esquerda)
PASTA_DADOS = "dados"
PASTA_FIGURAS = "figuras"
os.makedirs(PASTA_DADOS, exist_ok=True)
os.makedirs(PASTA_FIGURAS, exist_ok=True)

# Estilo dos gráficos: fonte grande e fundo claro, pensando em projetor
sns.set_theme(style="white", context="talk")
plt.rcParams.update({
    "font.size": 16, "axes.titlesize": 20, "axes.labelsize": 17,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14,
    "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
})

# Cores consistentes por grupo (paleta amigável a daltônicos)
COR_GRUPO = {
    "Primatas":      "#0072B2",
    "Carnívoros":    "#D55E00",
    "Ungulados e baleias": "#009E73",
    "Roedor":        "#CC79A7",
    "Ave (grupo externo)":   "#B8860B",
    "Peixe (grupo externo)": "#555555",
}"""))

# ---------------- 2. Dados
C.append(md("""## 2. Dados: baixando o gene CYTB do NCBI

Usamos o **Biopython (`Bio.Entrez`)** para baixar os registros do banco `nucleotide` do NCBI (GenBank).

**Como os accession numbers foram encontrados:** buscando no NCBI (`Entrez.esearch`) por organismo + gene CYTB + comprimento aproximado de 1.100–1.200 pb, e conferindo o título do registro (preferimos "complete cds"). Nenhum código foi inventado.
Para humano, gorila, camundongo e peixe-zebra não usamos um registro curto só do CYTB, e sim o **genoma mitocondrial de referência (RefSeq, `NC_...`)**, extraindo o trecho do CYTB pela anotação do próprio registro.

Cada sequência é **validada**: comprimento entre 1.080 e 1.200 bases e quase só letras A/C/G/T.
Se o NCBI estiver fora do ar (ou menos de 12 espécies forem baixadas), o notebook usa automaticamente um **conjunto menor embutido (Plano B)** e avisa claramente."""))

C.append(code("""# Tabela de espécies: (nome popular, nome científico, grupo, accession no NCBI)
ESPECIES = [
    ("Humano",        "Homo sapiens",               "Primatas",              "NC_012920.1"),   # genoma mito de referência
    ("Chimpanzé",     "Pan troglodytes",            "Primatas",              "X93338.1"),
    ("Gorila",        "Gorilla gorilla gorilla",    "Primatas",              "NC_011120.1"),   # genoma mito de referência
    ("Orangotango",   "Pongo abelii",               "Primatas",              "U38274.1"),
    ("Macaco-rhesus", "Macaca mulatta",             "Primatas",              "U38272.1"),
    ("Cachorro",      "Canis lupus familiaris",     "Carnívoros",            "EU352854.1"),
    ("Lobo",          "Canis lupus",                "Carnívoros",            "AY928668.1"),
    ("Gato",          "Felis catus",                "Carnívoros",            "KU253483.1"),
    ("Leão",          "Panthera leo",               "Carnívoros",            "AF384818.1"),
    ("Vaca",          "Bos taurus",                 "Ungulados e baleias",   "PP997875.1"),
    ("Hipopótamo",    "Hippopotamus amphibius",     "Ungulados e baleias",   "U07565.1"),
    ("Baleia-azul",   "Balaenoptera musculus",      "Ungulados e baleias",   "KJ586852.1"),
    ("Camundongo",    "Mus musculus",               "Roedor",                "NC_005089.1"),   # genoma mito de referência
    ("Galinha",       "Gallus gallus",              "Ave (grupo externo)",   "AB044986.1"),
    ("Peixe-zebra",   "Danio rerio",                "Peixe (grupo externo)", "NC_002333.2"),   # genoma mito de referência
]
GRUPO_DE = {nome: grupo for nome, _, grupo, _ in ESPECIES}
CIENT_DE = {nome: cient for nome, cient, _, _ in ESPECIES}

COMPRIMENTO_MIN, COMPRIMENTO_MAX = 1080, 1200   # "próximo de 1.140 bases"
MIN_ESPECIES = 12      # abaixo disso, usamos o Plano B
FORCAR_PLANO_B = False # coloque True para testar o Plano B mesmo com internet


def baixar_cytb(accession, tentativas=3):
    \"\"\"Baixa um registro do NCBI e devolve a sequência do CYTB (str).
    - Se o registro tiver uma feature CDS de CYTB (caso dos genomas mitocondriais), extrai só esse trecho.
    - Caso contrário (registros que já são só o gene), usa a sequência inteira.\"\"\"
    ultimo_erro = None
    for t in range(tentativas):
        try:
            h = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            reg = SeqIO.read(h, "genbank")
            h.close()
            time.sleep(0.4)                       # o NCBI limita a ~3 requisições por segundo
            for f in reg.features:
                if f.type == "CDS":
                    txt = (f.qualifiers.get("gene", [""])[0] + " " + f.qualifiers.get("product", [""])[0]).lower()
                    if any(k in txt for k in ("cytb", "cyt b", "cytochrome b", "cob")):
                        return str(f.extract(reg.seq)).upper()
            return str(reg.seq).upper()
        except Exception as e:                    # rede instável, limite de uso, etc.
            ultimo_erro = e
            time.sleep(1.5 * (t + 1))
    raise RuntimeError(f"falha ao baixar {accession}: {ultimo_erro}")


def validar(seq):
    \"\"\"Valida comprimento (~1.140 pb) e composição (A/C/G/T). Devolve (ok, motivo).\"\"\"
    if not (COMPRIMENTO_MIN <= len(seq) <= COMPRIMENTO_MAX):
        return False, f"comprimento {len(seq)} fora de {COMPRIMENTO_MIN}-{COMPRIMENTO_MAX}"
    frac_acgt = sum(seq.count(b) for b in "ACGT") / len(seq)
    if frac_acgt < 0.99:
        return False, f"só {frac_acgt:.1%} das letras são A/C/G/T"
    return True, "OK"


# ---------- Download com validação ----------
baixadas, status = {}, {}
if not FORCAR_PLANO_B:
    for nome, cient, grupo, acc in ESPECIES:
        try:
            s = baixar_cytb(acc)
            ok, motivo = validar(s)
            status[nome] = motivo
            if ok:
                baixadas[nome] = s
        except Exception as e:
            status[nome] = f"ERRO: {e}"[:80]
        print(f"{nome:14s} {acc:12s} -> {status[nome]}")
print(f"\\n{len(baixadas)} de {len(ESPECIES)} espécies baixadas e validadas.")"""))

C.append(md("""### Plano B: sequências embutidas
Esta célula só é usada se o NCBI estiver inacessível. Contém **9 espécies** (humano, chimpanzé, gorila, cachorro, gato, vaca, baleia-azul, galinha e peixe-zebra), copiadas das mesmas sequências do NCBI baixadas durante o desenvolvimento do notebook (accession no cabeçalho de cada uma)."""))

C.append(code('''# Plano B: sequências embutidas (FASTA). Só será usado se o download falhar.
FASTA_PLANO_B = """
''' + PLANO_B + '''
"""'''))

C.append(code("""# ---------- Decide entre dados do NCBI e Plano B ----------
if len(baixadas) >= MIN_ESPECIES and not FORCAR_PLANO_B:
    FONTE = "NCBI"
    seqs = baixadas
    ACC_DE = {nome: acc for nome, _, _, acc in ESPECIES}
    ordem = [n for n, *_ in ESPECIES if n in seqs]
else:
    FONTE = "PLANO B (sequências embutidas)"
    print("=" * 78)
    print("ATENÇÃO: usando o PLANO B — conjunto reduzido de sequências embutidas no notebook.")
    print("Motivo: NCBI inacessível, download incompleto ou FORCAR_PLANO_B=True.")
    print("=" * 78)
    seqs, ACC_DE = {}, {}
    for r in SeqIO.parse(io.StringIO(FASTA_PLANO_B.strip()), "fasta"):
        nome = r.description.split(" | ")[0].split(" ", 1)[1]
        seqs[nome] = str(r.seq).upper()
        ACC_DE[nome] = r.id
    ordem = [n for n, *_ in ESPECIES if n in seqs]

# ---------- Salva em FASTA ----------
registros = [SeqRecord(Seq(seqs[n]), id=ACC_DE[n], description=f"{n} | {CIENT_DE[n]} | CYTB") for n in ordem]
caminho_fasta = os.path.join(PASTA_DADOS, "cytb_especies.fasta")
SeqIO.write(registros, caminho_fasta, "fasta")

# ---------- Tabela (espécie, accession, comprimento) ----------
tabela = pd.DataFrame({
    "Espécie": ordem,
    "Nome científico": [CIENT_DE[n] for n in ordem],
    "Grupo": [GRUPO_DE[n] for n in ordem],
    "Accession": [ACC_DE[n] for n in ordem],
    "Comprimento (pb)": [len(seqs[n]) for n in ordem],
    "Validação": [validar(seqs[n])[1] for n in ordem],
})
print(f"Fonte dos dados: {FONTE}")
print(f"FASTA salvo em: {caminho_fasta}   ({len(ordem)} sequências)")
display(tabela)"""))

# ---------------- 3. Algoritmo
C.append(md("""## 3. O algoritmo: alinhamento global de Needleman-Wunsch (do zero)

**A ideia.** Montamos uma tabela `F` com `(n+1) × (m+1)` células, onde a célula `F[i][j]` guarda **a melhor pontuação possível ao alinhar os primeiros `i` caracteres da sequência A com os primeiros `j` da sequência B**.
Cada célula só pode vir de 3 vizinhas:

| Movimento | Significado | Pontuação |
|---|---|---|
| ↖ diagonal | alinha `A[i]` com `B[j]` | `F[i-1][j-1] + (match ou mismatch)` |
| ↑ cima | `A[i]` fica contra um buraco | `F[i-1][j] + gap` |
| ← esquerda | `B[j]` fica contra um buraco | `F[i][j-1] + gap` |

`F[i][j]` = o **maior** dos três. A primeira linha e a primeira coluna são só buracos (`i·gap` e `j·gap`).
No fim, o **traceback** parte do canto inferior direito e volta seguindo de onde cada valor veio, reconstruindo o alinhamento.

Os parâmetros de pontuação (`match`, `mismatch`, `gap`) podem ser alterados."""))

C.append(code("""# ===== Needleman-Wunsch, versão didática (Python puro) =====
# Pontuação padrão (parametrizável em todas as funções)
MATCH, MISMATCH, GAP = 1, -1, -2


def nw_matriz_python(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP):
    \"\"\"Preenche a matriz de pontuação F com dois laços for (O(n*m)).\"\"\"
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):              # borda: alinhar prefixo de A contra só buracos
        F[i][0] = i * gap
    for j in range(1, m + 1):              # borda: alinhar prefixo de B contra só buracos
        F[0][j] = j * gap
    for i in range(1, n + 1):
        ai, Fi, Fa = a[i - 1], F[i], F[i - 1]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            # melhor entre: diagonal, vir de cima, vir da esquerda
            Fi[j] = max(Fa[j - 1] + s, Fa[j] + gap, Fi[j - 1] + gap)
    return np.array(F, dtype=np.int32)


def nw_traceback(F, a, b, match=MATCH, mismatch=MISMATCH, gap=GAP):
    \"\"\"Reconstrói um alinhamento ótimo, andando do canto (n,m) até (0,0).
    Em caso de empate, prefere a diagonal. Devolve (a_alinhada, b_alinhada, caminho).\"\"\"
    i, j = len(a), len(b)
    al_a, al_b, caminho = [], [], [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i, j] == F[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            al_a.append(a[i - 1]); al_b.append(b[j - 1]); i -= 1; j -= 1      # diagonal
        elif i > 0 and F[i, j] == F[i - 1, j] + gap:
            al_a.append(a[i - 1]); al_b.append("-"); i -= 1                   # buraco em B
        else:
            al_a.append("-"); al_b.append(b[j - 1]); j -= 1                   # buraco em A
        caminho.append((i, j))
    return "".join(reversed(al_a)), "".join(reversed(al_b)), caminho[::-1]


def identidade(al_a, al_b):
    \"\"\"Identidade percentual = colunas iguais / total de colunas do alinhamento (buracos contam como diferença).\"\"\"
    iguais = sum(1 for x, y in zip(al_a, al_b) if x == y and x != "-")
    return 100.0 * iguais / len(al_a)


def mostrar_alinhamento(al_a, al_b, largura=60):
    \"\"\"Imprime o alinhamento em blocos, com '|' nas colunas iguais.\"\"\"
    marca = "".join("|" if x == y and x != "-" else " " for x, y in zip(al_a, al_b))
    for k in range(0, len(al_a), largura):
        print(al_a[k:k + largura]); print(marca[k:k + largura]); print(al_b[k:k + largura]); print()


def alinhar(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP):
    \"\"\"Atalho: matriz + traceback + identidade (versão Python puro).\"\"\"
    F = nw_matriz_python(a, b, match, mismatch, gap)
    al_a, al_b, caminho = nw_traceback(F, a, b, match, mismatch, gap)
    return F, al_a, al_b, caminho, identidade(al_a, al_b)"""))

C.append(md("""### 3.1 Teste com um exemplo curto
Primeiro as frases do exemplo da introdução (letras no lugar de bases) e depois dois trechos curtos de DNA, com a matriz de pontuação visualizada. O caminho do traceback aparece destacado."""))

C.append(code("""# Teste 1: as frases do exemplo (o algoritmo funciona com qualquer texto!)
F, al_a, al_b, _, ident = alinhar("O GATO COMEU O PEIXE", "O GATO COMEU PEIXE")
mostrar_alinhamento(al_a, al_b)
print(f"Pontuação final: {F[-1, -1]}   |   identidade: {ident:.1f}%")"""))

C.append(code("""def plotar_matriz(F, a, b, caminho, arquivo=None):
    \"\"\"Desenha a matriz de pontuação com os valores e o caminho do traceback destacado.\"\"\"
    fig, ax = plt.subplots(figsize=(9.5, 8))
    sns.heatmap(F, annot=True, fmt="d", cmap="Blues", cbar_kws={"label": "Pontuação"},
                xticklabels=["-"] + list(b), yticklabels=["-"] + list(a),
                linewidths=0.5, linecolor="white", annot_kws={"size": 13}, ax=ax)
    ys, xs = zip(*caminho)
    ax.plot(np.array(xs) + 0.5, np.array(ys) + 0.5, color="#D55E00", lw=4, marker="o", ms=9, label="Caminho do traceback")
    ax.set_xlabel("Sequência B"); ax.set_ylabel("Sequência A")
    ax.set_title("Matriz de pontuação de Needleman-Wunsch", pad=40)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), frameon=True)
    ax.xaxis.tick_top(); ax.xaxis.set_label_position("top")
    if arquivo:
        fig.savefig(os.path.join(PASTA_FIGURAS, arquivo))
    plt.show()


# Teste 2: dois trechos curtos de DNA (match=+1, mismatch=-1, gap=-2)
seq_a, seq_b = "GATTACA", "GCATGCA"
F_ex, al_a, al_b, caminho, ident = alinhar(seq_a, seq_b)
mostrar_alinhamento(al_a, al_b)
print(f"Pontuação final: {F_ex[-1, -1]}   |   identidade: {ident:.1f}%")
plotar_matriz(F_ex, seq_a, seq_b, caminho, arquivo="01_matriz_needleman_wunsch.png")"""))

C.append(md("""### 3.2 Versão otimizada com NumPy
O algoritmo em Python puro faz ~1,3 milhão de cálculos por par de genes (1.140 × 1.140), e precisamos de **105 pares** (15 espécies). Dá para otimizar:

* **Vetorização por linhas.** Em cada linha `i` calculamos todas as colunas de uma vez com NumPy: a diagonal e o "vir de cima" dependem só da linha anterior, então viram operações de vetor.
* **O único problema** é o "vir da esquerda" (`F[i][j-1] + gap`), que depende da célula ao lado. Como o gap é linear, isso equivale a `F[i][j] = j·gap + máximo acumulado de (H[k] − k·gap)` para `k ≤ j`, onde `H` é o melhor entre diagonal e cima. Isso é um `np.maximum.accumulate`, que o NumPy executa em C.
* **Pré-computamos** o vetor de match/mismatch de cada letra (A, C, G, T) contra a sequência B, em vez de comparar caracteres um a um.

O resultado é **exatamente a mesma matriz** (verificamos abaixo), só que bem mais rápido. O traceback continua o mesmo, pois custa apenas O(n+m)."""))

C.append(code("""def nw_matriz_numpy(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP):
    \"\"\"Mesma matriz F do Needleman-Wunsch, mas calculada linha a linha com NumPy.\"\"\"
    n, m = len(a), len(b)
    b_arr = np.frombuffer(b.encode(), dtype=np.uint8)
    # vetor de pontuação de cada letra contra todas as posições de B (calculado uma única vez)
    pont = {c: np.where(b_arr == ord(c), match, mismatch).astype(np.int32) for c in set(a)}
    kgap = np.arange(m + 1, dtype=np.int32) * gap          # [0, gap, 2*gap, ...]
    F = np.empty((n + 1, m + 1), dtype=np.int32)
    F[0] = kgap                                            # primeira linha: só buracos
    E = np.empty(m + 1, dtype=np.int32)
    for i in range(1, n + 1):
        E[0] = i * gap                                     # primeira coluna: só buracos
        E[1:] = np.maximum(F[i - 1, :-1] + pont[a[i - 1]], # diagonal (match/mismatch)
                           F[i - 1, 1:] + gap)             # vindo de cima
        # "vindo da esquerda" resolvido com máximo acumulado (substitui o laço interno)
        F[i] = np.maximum.accumulate(E - kgap) + kgap
    return F


def alinhar_rapido(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP):
    \"\"\"Mesmo que alinhar(), usando a matriz vetorizada.\"\"\"
    F = nw_matriz_numpy(a, b, match, mismatch, gap)
    al_a, al_b, caminho = nw_traceback(F, a, b, match, mismatch, gap)
    return F, al_a, al_b, caminho, identidade(al_a, al_b)


# Verificação: as duas versões produzem a MESMA matriz (em exemplos curtos e em genes reais)
assert np.array_equal(nw_matriz_python(seq_a, seq_b), nw_matriz_numpy(seq_a, seq_b))
x, y = seqs[ordem[0]][:300], seqs[ordem[-1]][:300]
assert np.array_equal(nw_matriz_python(x, y), nw_matriz_numpy(x, y))
x, y = seqs[ordem[0]][:300], seqs[ordem[1]][:290]
assert np.array_equal(nw_matriz_python(x, y, 2, -3, -4), nw_matriz_numpy(x, y, 2, -3, -4))
print("OK: versão NumPy idêntica à versão em Python puro (inclusive com outros parâmetros de pontuação).")"""))

C.append(md("""### 3.3 Validação contra o `PairwiseAligner` do Biopython
Configuramos o Biopython com a **mesma pontuação** (global, match +1, mismatch −1, gap −2 linear) e comparamos, em pares de genes completos:
* a **pontuação ótima** (deve ser idêntica; é o que prova que o algoritmo está correto);
* a **identidade percentual** (pode diferir um pouco, porque podem existir vários alinhamentos igualmente ótimos e cada programa escolhe um deles no desempate)."""))

C.append(code("""# Aligner do Biopython com a mesma pontuação do nosso algoritmo
aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score, aligner.mismatch_score = MATCH, MISMATCH
aligner.open_gap_score, aligner.extend_gap_score = GAP, GAP     # gap linear (igual ao nosso)

pares_teste = [("Humano", "Chimpanzé"), ("Humano", "Gorila"), ("Cachorro", "Lobo"),
               ("Gato", "Leão"), ("Baleia-azul", "Hipopótamo"), ("Humano", "Galinha"), ("Humano", "Peixe-zebra")]
pares_teste = [(p, q) for p, q in pares_teste if p in seqs and q in seqs]

linhas = []
for p, q in pares_teste:
    F, al_a, al_b, _, ident = alinhar_rapido(seqs[p], seqs[q])
    meu_score = int(F[-1, -1])
    bio_score = aligner.score(seqs[p], seqs[q])
    aln = aligner.align(seqs[p], seqs[q])[0]
    bio_ident = 100.0 * aln.counts().identities / aln.length
    linhas.append({"Par": f"{p} × {q}", "Score (meu)": meu_score, "Score (Biopython)": int(bio_score),
                   "Scores iguais?": meu_score == int(bio_score),
                   "Identidade % (meu)": round(ident, 2), "Identidade % (Biopython)": round(bio_ident, 2)})
comparacao = pd.DataFrame(linhas)
display(comparacao)
assert comparacao["Scores iguais?"].all(), "Pontuações divergem do Biopython!"
print("Validado: a pontuação ótima do meu algoritmo coincide com a do Biopython em todos os pares testados.")

# Mostra um alinhamento real (primeiras 120 colunas) para dar uma ideia de como fica
p, q = pares_teste[0]
_, al_a, al_b, _, ident = alinhar_rapido(seqs[p], seqs[q])
print(f"\\nTrecho do alinhamento {p} × {q} (identidade total {ident:.1f}%):\\n")
mostrar_alinhamento(al_a[:120], al_b[:120])"""))

# ---------------- 4. Análise
C.append(md("""## 4. Análise: matriz de distância, heatmap e árvore

Para cada par de espécies:
1. alinhamos os genes (Needleman-Wunsch, versão NumPy);
2. calculamos a **identidade** (fração de colunas iguais no alinhamento);
3. definimos a **distância = 1 − identidade**. Quanto menor, mais parecidos os DNAs.

Só calculamos o triângulo superior (a distância de A a B é igual à de B a A) e espelhamos."""))

C.append(code("""def matriz_distancia(seqs, nomes, alinhador=alinhar_rapido):
    \"\"\"Distância = 1 - identidade para todos os pares (só triângulo superior; a matriz é simétrica).\"\"\"
    k = len(nomes)
    D = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            *_, ident = alinhador(seqs[nomes[i]], seqs[nomes[j]])
            D[i, j] = D[j, i] = 1.0 - ident / 100.0
    return D


# Estimativa do ganho da otimização: tempo de UM par nas duas versões
a, b = seqs[ordem[0]], seqs[ordem[1]]
t0 = time.perf_counter(); alinhar(a, b); t_py = time.perf_counter() - t0
t0 = time.perf_counter(); alinhar_rapido(a, b); t_np = time.perf_counter() - t0
n_pares = len(ordem) * (len(ordem) - 1) // 2
print(f"Um par (~1.140 × 1.140 bases):  Python puro = {t_py:.2f} s   |   NumPy = {t_np:.3f} s   (≈ {t_py / t_np:.0f}x mais rápido)")
print(f"Para {n_pares} pares: Python puro ≈ {t_py * n_pares / 60:.1f} min   |   NumPy ≈ {t_np * n_pares:.0f} s")

t0 = time.perf_counter()
D = matriz_distancia(seqs, ordem)
print(f"\\nMatriz de distância calculada em {time.perf_counter() - t0:.1f} s (versão NumPy).")

dist_df = pd.DataFrame(D, index=ordem, columns=ordem)
dist_df.to_csv(os.path.join(PASTA_DADOS, "matriz_distancia.csv"), float_format="%.4f")
display((dist_df * 100).round(1).style.background_gradient(cmap="viridis_r").format("{:.1f}"))"""))

C.append(code("""# ---------- UPGMA: linkage "average" do scipy ----------
# UPGMA = agrupar sempre os dois grupos mais próximos, usando a MÉDIA das distâncias entre os membros.
Z = linkage(squareform(D, checks=False), method="average")
folhas = dendrogram(Z, no_plot=True)["leaves"]          # ordem das espécies na árvore
ordem_arvore = [ordem[i] for i in folhas]

def colorir_rotulos(ax, eixo="y"):
    \"\"\"Pinta o texto de cada espécie com a cor do seu grupo.\"\"\"
    rotulos = ax.get_yticklabels() if eixo == "y" else ax.get_xticklabels()
    for r in rotulos:
        r.set_color(COR_GRUPO[GRUPO_DE[r.get_text()]]); r.set_fontweight("bold")

def legenda_grupos(ax, loc="lower right", **kw):
    grupos = [g for g in COR_GRUPO if g in {GRUPO_DE[n] for n in ordem}]
    ax.legend(handles=[Patch(color=COR_GRUPO[g], label=g) for g in grupos], loc=loc, frameon=True, title="Grupo", **kw)

# ---------- (a) Heatmap (espécies na ordem da árvore, para os blocos aparecerem) ----------
idx = [ordem.index(n) for n in ordem_arvore]
D_ord = D[np.ix_(idx, idx)] * 100
fig, ax = plt.subplots(figsize=(15, 12.5))
sns.heatmap(D_ord, annot=True, fmt=".0f", cmap="viridis_r", vmin=0, square=True,
            xticklabels=ordem_arvore, yticklabels=ordem_arvore, linewidths=0.4, linecolor="white",
            annot_kws={"size": 12}, cbar_kws={"label": "Diferença entre os genes (%) = (1 − identidade) × 100", "shrink": 0.8}, ax=ax)
ax.set_title("Quão diferentes são os genes CYTB? (amarelo = mais parecidos)")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
colorir_rotulos(ax, "x"); colorir_rotulos(ax, "y")
legenda_grupos(ax, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)   # legenda fora da matriz
fig.savefig(os.path.join(PASTA_FIGURAS, "02_heatmap_distancias.png"))
plt.show()"""))

C.append(code("""# ---------- (b) Árvore (dendrograma UPGMA) ----------
fig, ax = plt.subplots(figsize=(13, 9))
dendrogram(Z, labels=ordem, orientation="right", ax=ax, link_color_func=lambda k: "#333333", leaf_font_size=17)
ax.set_title("Árvore de parentesco (UPGMA) a partir do gene CYTB")
ax.set_xlabel("Distância média (1 − identidade): junção mais à esquerda = parentesco mais próximo")
colorir_rotulos(ax, "y")
for lado in ("top", "right"):
    ax.spines[lado].set_visible(False)
ax.grid(axis="x", alpha=0.3)
legenda_grupos(ax, loc="upper right")
fig.savefig(os.path.join(PASTA_FIGURAS, "03_arvore_upgma.png"))
plt.show()"""))

C.append(md("""### O que esperávamos ver
Pela biologia já conhecida, esperamos que:
* **primatas** fiquem juntos (e, dentro deles, humano–chimpanzé–gorila mais próximos entre si);
* **felinos** (gato e leão) fiquem juntos e **canídeos** (cachorro e lobo) também, formando um grupo de carnívoros;
* a **baleia** fique perto do **hipopótamo** (os hipopótamos são os parentes vivos mais próximos das baleias!) e ambos perto da vaca;
* **galinha** e **peixe-zebra** (grupo externo) fiquem bem afastados dos mamíferos, com o peixe como o mais distante.

Vamos conferir isso de forma automática, perguntando para a árvore se cada grupo forma um ramo próprio (um *clado*)."""))

C.append(code("""def e_clado(Z, nomes, grupo):
    \"\"\"True se as espécies do grupo formam um ramo exclusivo da árvore (e nenhuma outra espécie está nele).\"\"\"
    _, nos = to_tree(Z, rd=True)
    alvo = {nomes.index(n) for n in grupo}
    return any(set(no.pre_order()) == alvo for no in nos)

def vizinho_mais_proximo(nome):
    i = ordem.index(nome)
    d = D[i].copy(); d[i] = np.inf
    return ordem[int(np.argmin(d))], d.min()

testes = {
    "Primatas formam um ramo":                 [n for n in ordem if GRUPO_DE[n] == "Primatas"],
    "Humano + chimpanzé + gorila juntos":      ["Humano", "Chimpanzé", "Gorila"],
    "Gato + leão juntos (felinos)":            ["Gato", "Leão"],
    "Cachorro + lobo juntos (canídeos)":       ["Cachorro", "Lobo"],
    "Carnívoros formam um ramo":               [n for n in ordem if GRUPO_DE[n] == "Carnívoros"],
    "Baleia + hipopótamo juntos":              ["Baleia-azul", "Hipopótamo"],
    "Baleia + hipopótamo + vaca juntos":       ["Baleia-azul", "Hipopótamo", "Vaca"],
    "Todos os mamíferos formam um ramo (grupo externo fora)": [n for n in ordem if GRUPO_DE[n] not in ("Ave (grupo externo)", "Peixe (grupo externo)")],
}
print("Conferência da árvore com o esperado:\\n")
for desc, grupo in testes.items():
    if len(grupo) >= 2 and all(n in ordem for n in grupo):
        print(f"  {'✔' if e_clado(Z, ordem, grupo) else '✘'}  {desc}")

print("\\nVizinho mais próximo de cada espécie (menor distância):")
for n in ordem:
    v, d = vizinho_mais_proximo(n)
    print(f"  {n:14s} → {v:14s} (diferença de {d * 100:.1f}%)")"""))

C.append(md("""__INTERPRETACAO__"""))

# ---------------- 5. Complexidade
C.append(md("""## 5. Complexidade: quanto tempo o algoritmo leva?

O preenchimento da matriz calcula **uma célula para cada par (i, j)**, e cada célula custa um número fixo de operações. Logo o tempo é proporcional a **n × m**, ou seja, **O(n·m)**. Com sequências do mesmo tamanho `n`, isso vira **O(n²)**: dobrar o tamanho deve **quadruplicar** o tempo.
(O traceback custa só O(n+m), desprezível perto do preenchimento; a memória também é O(n·m), por guardarmos a matriz inteira.)

Vamos medir: alinhamos prefixos de tamanhos crescentes (humano × vaca) e comparamos Python puro, NumPy e o Biopython (em C)."""))

C.append(code("""tamanhos = [100, 200, 400, 800, 1100]
tamanhos = [t for t in tamanhos if t <= min(len(s) for s in seqs.values())]
p, q = ("Humano", "Vaca") if {"Humano", "Vaca"} <= set(seqs) else (ordem[0], ordem[-1])

def melhor_tempo(f, repeticoes):
    \"\"\"Menor tempo entre várias repetições (reduz o ruído de outros processos).\"\"\"
    ts = []
    for _ in range(repeticoes):
        t0 = time.perf_counter(); f(); ts.append(time.perf_counter() - t0)
    return min(ts)

res = []
for n in tamanhos:
    a, b = seqs[p][:n], seqs[q][:n]
    res.append({
        "n (bases)": n,
        "células (n·m)": n * n,
        "Python puro (s)": melhor_tempo(lambda: nw_matriz_python(a, b), 3 if n <= 400 else 2),
        "NumPy (s)":       melhor_tempo(lambda: nw_matriz_numpy(a, b), 5),
        "Biopython (s)":   melhor_tempo(lambda: aligner.score(a, b), 5),
    })
tempos = pd.DataFrame(res)
tempos["Python puro: µs por célula"] = tempos["Python puro (s)"] / tempos["células (n·m)"] * 1e6
display(tempos.style.format({"Python puro (s)": "{:.4f}", "NumPy (s)": "{:.4f}", "Biopython (s)": "{:.5f}",
                             "Python puro: µs por célula": "{:.3f}", "células (n·m)": "{:,}"}))
tempos.to_csv(os.path.join(PASTA_DADOS, "tempos_complexidade.csv"), index=False)

# Expoente empírico: inclinação da reta no gráfico log-log (tempo ~ n^expoente)
for col in ["Python puro (s)", "NumPy (s)", "Biopython (s)"]:
    expo = np.polyfit(np.log(tempos["n (bases)"]), np.log(tempos[col]), 1)[0]
    print(f"Expoente medido para {col:16s}: tempo ∝ n^{expo:.2f}")"""))

C.append(code("""fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 6.8))
estilos = [("Python puro (s)", "#D55E00", "o", "Python puro (meu algoritmo)"),
           ("NumPy (s)",      "#0072B2", "s", "NumPy (meu algoritmo otimizado)"),
           ("Biopython (s)",  "#009E73", "^", "Biopython (C)")]

# Painel 1: escala linear, o "estouro" do Python puro
for col, cor, marc, rot in estilos[:2]:
    ax1.plot(tempos["n (bases)"], tempos[col], marker=marc, color=cor, lw=3, ms=10, label=rot)
ax1.set_title("Tempo × tamanho (escala linear)")
ax1.set_xlabel("Tamanho de cada sequência (bases)"); ax1.set_ylabel("Tempo (segundos)")
ax1.legend(); ax1.grid(alpha=0.3)

# Painel 2: log-log; uma reta de inclinação 2 significa O(n²) = O(n·m)
for col, cor, marc, rot in estilos:
    ax2.plot(tempos["n (bases)"], tempos[col], marker=marc, color=cor, lw=3, ms=10, label=rot)
n0, t0 = tempos["n (bases)"].iloc[0], tempos["Python puro (s)"].iloc[0]
ns = np.array(tempos["n (bases)"])
ax2.plot(ns, t0 * (ns / n0) ** 2, "k--", lw=2, label="Referência: crescimento O(n²)")
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xticks(tamanhos); ax2.set_xticklabels(tamanhos); ax2.minorticks_off()
ax2.set_title("Mesmos dados em escala log-log")
ax2.set_xlabel("Tamanho de cada sequência (bases)"); ax2.set_ylabel("Tempo (segundos)")
ax2.legend(fontsize=12); ax2.grid(alpha=0.3, which="both")
fig.suptitle("Complexidade do alinhamento global: o tempo cresce como n·m", fontsize=22, y=1.02)
fig.savefig(os.path.join(PASTA_FIGURAS, "04_complexidade_tempo.png"))
plt.show()"""))

C.append(md("""__COMPLEXIDADE_TXT__"""))

# ---------------- 6. Conclusões
C.append(md("""__CONCLUSOES__"""))

C.append(code("""# Resumo dos arquivos gerados
print("Arquivos gerados:")
for pasta in (PASTA_DADOS, PASTA_FIGURAS):
    for f in sorted(os.listdir(pasta)):
        print(f"  {pasta}/{f}  ({os.path.getsize(os.path.join(pasta, f)) / 1024:.0f} KB)")
print(f"\\nFonte dos dados desta execução: {FONTE}")"""))

# ---- textos finais (preenchidos no segundo passo, se existirem)
TXT = {}
for k in ("INTERPRETACAO", "COMPLEXIDADE_TXT", "CONCLUSOES"):
    p = os.path.join(AQUI, f"txt_{k}.md")
    TXT[k] = open(p, encoding="utf-8").read() if os.path.exists(p) else f"*(texto {k} pendente)*"
for c in C:
    for k, v in TXT.items():
        c.source = c.source.replace(f"__{k}__", v)

nb = new_notebook(cells=C)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
nb.metadata["colab"] = {"provenance": []}
saida = os.path.join(AQUI, "projeto_final_semsaida.ipynb")
nbformat.write(nb, saida)
print("gerado:", saida)
