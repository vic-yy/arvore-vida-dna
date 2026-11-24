# Roteiro de apresentação (15 min): "Quem é parente de quem? A árvore da vida pelo DNA"

Público misto, com pessoas leigas. A fala é sugerida, então adapte ao seu jeito. Os tempos somam **15 min**; o bloco 4 é o principal.
Figuras de apoio estão na pasta `figuras/`.

| Bloco | Tema | Tempo | Slides |
|---|---|---|---|
| 1 | Pergunta e DNA como texto | 2 min | 1–2 |
| 2 | Os dados | 2 min | 3 |
| 3 | O algoritmo (alinhar sequências) | 4 min | 4–6 |
| 4 | O resultado: heatmap e árvore | 4 min | 7–9 |
| 5 | Custo computacional | 2 min | 10 |
| 6 | Conclusões e limitações | 1 min | 11 |

---

## Bloco 1 — A pergunta (0:00–2:00)

**Slide 1 – Título + pergunta.** Mostrar: título, e uma foto de humano, chimpanzé, cachorro e baleia. Pergunta grande: *"Dá para descobrir quem é parente de quem só olhando o DNA?"*

**Slide 2 – DNA é um texto.** Mostrar: um trecho de DNA (`ATGACCCCAATACGC…`) ao lado de uma frase em português. Legenda: *4 letras: A, C, G, T*.

> **Fala:** "Pensem num primo distante: vocês dividem avós, mas já se parecem menos do que dois irmãos. Com as espécies é igual: quanto mais antigo o ancestral em comum, mais diferenças se acumulam. O DNA é um texto escrito com só quatro letras, A, C, G e T. Se eu comparar o 'texto' de dois animais e contar quantas letras diferem, tenho uma régua de parentesco. Hoje vou usar um gene de cerca de 1.140 letras, o citocromo b, que fica nas mitocôndrias, as usinas de energia da célula."

## Bloco 2 — Os dados (2:00–4:00)

**Slide 3 – 15 espécies, 1 gene.** Mostrar: tabela (espécie, accession) do notebook ou os 15 animais coloridos por grupo (primatas, carnívoros, ungulados e baleias, roedor, ave, peixe). Destacar: *galinha e peixe-zebra = grupo externo*.

> **Fala:** "Baixei esse gene de 15 espécies direto do GenBank, o grande banco público de DNA do NCBI. Escolhi primatas, cães e gatos, vaca, hipopótamo e baleia, um camundongo, e dois 'penetras' de propósito: uma galinha e um peixe. Eles servem de controle: se meu método funcionar, esses dois têm que ficar bem longe dos mamíferos. Verifiquei que todas as sequências têm perto de 1.140 letras. Todos eles têm um gene com o mesmo papel, então dá para comparar um com o outro."

## Bloco 3 — O algoritmo (4:00–8:00)

**Slide 4 – Analogia das frases.** Mostrar:
```
O GATO COMEU O PEIXE
O GATO COMEU - PEIXE
```
**Slide 5 – Matriz de pontuação** (`01_matriz_needleman_wunsch.png`). Destacar o caminho laranja na diagonal.
**Slide 6 – Validação.** Mostrar a tabela de comparação com o Biopython: *scores iguais em todos os pares*.

> **Fala (analogia):** "Como comparar dois textos que não são idênticos? Olhem essas duas frases. Se eu simplesmente comparar letra por letra, a partir da palavra que falta tudo fica deslocado e parece que as frases são totalmente diferentes. **Alinhar** é deslizar as frases e deixar um buraco onde falta algo, para que as partes parecidas fiquem uma embaixo da outra. Aí conto: 18 colunas iguais de 20, ou seja, 90% de identidade. No DNA é a mesma coisa: um buraco é uma letra que foi inserida ou perdida na evolução, e uma coluna diferente é uma mutação."
>
> "Como o computador acha o *melhor* alinhamento, entre bilhões de possibilidades? Com o algoritmo de Needleman-Wunsch, que eu programei do zero. Ele monta uma tabela como esta: cada casinha guarda a melhor pontuação possível até aquele ponto. Letra igual dá ponto, letra diferente ou buraco tira ponto. Depois, é só seguir o caminho de volta, a linha laranja, e o alinhamento aparece."
>
> "E como saber se meu código está certo? Comparei com o Biopython, uma biblioteca profissional. Nos 7 pares que testei, a pontuação foi idêntica. Também criei uma versão otimizada com NumPy, que dá exatamente o mesmo resultado, mas bem mais rápido: um par de genes em cerca de 0,01 segundo."

## Bloco 4 — O resultado (8:00–12:00)

**Slide 7 – Matriz/heatmap** (`02_heatmap_distancias.png`). Destacar os blocos amarelo-esverdeados: primatas, cachorro-lobo, gato-leão.
**Slide 8 – Árvore** (`03_arvore_upgma.png`). **Slide principal.** Destacar primatas, carnívoros e a galinha/peixe isolados.
**Slide 9 – Bateu ou não?** Mostrar a lista ✔/✘ do notebook.

> **Fala:** "Para cada par de espécies, medi a distância como 1 menos a identidade: quanto menor, mais parecidos. Dá uma tabela 15 por 15, que eu pintei: quanto mais amarelo, mais parecidos os genes. Reparem nos quadrados claros ao longo da diagonal: são os grupos de parentes."
>
> "Agora a árvore. O método, chamado UPGMA, funciona assim: junta os dois mais parecidos, trata os dois como um só, e repete. O que ela mostrou? **Os primatas ficaram juntos**, com humano e chimpanzé mais próximos, depois gorila, orangotango e macaco. **Cachorro e lobo são quase idênticos**, só 0,4% de diferença, porque o cachorro é um lobo domesticado. **Gato e leão** formam outro par, e os dois se juntam como carnívoros. E a galinha e o peixe, que eram meus penetras, ficaram lá no fundo, com o peixe como o mais distante. Ou seja: só comparando letras, o computador redescobriu a árvore da vida."
>
> "Mas nem tudo bateu, e isso é parte da graça. Esperava a baleia colada no hipopótamo, que é seu parente vivo mais próximo. A árvore juntou o hipopótamo primeiro à vaca, mas a diferença é de só 0,2 ponto, praticamente um empate. E o camundongo, que é parente de primatas, apareceu perto dos carnívoros, porque o DNA mitocondrial dos roedores muda muito rápido. Um gene só não é suficiente para os parentescos muito antigos."

## Bloco 5 — Custo computacional (12:00–14:00)

**Slide 10 – Tempo × tamanho** (`04_complexidade_tempo.png`). Gráfico da direita (log-log), com a linha tracejada O(n²).

> **Fala:** "Quanto custa fazer isso? A tabela tem uma casinha para cada par de letras, então dobrar o tamanho das sequências **quadruplica** o tempo. Medi isso com sequências de 100 a 1.100 letras: o meu código em Python puro acompanha direitinho a linha tracejada do crescimento ao quadrado. A versão com NumPy é cerca de 20 vezes mais rápida, e a do Biopython, escrita em C, é a mais rápida de todas. Para 15 espécies isso é rapidinho, mas para genomas gigantes esse algoritmo seria inviável, e por isso os pesquisadores usam atalhos."

## Bloco 6 — Conclusões e limitações (14:00–15:00)

**Slide 11 – Conclusões.** Três tópicos: ✔ DNA funciona como régua de parentesco; ✔ algoritmo próprio validado contra o Biopython; ⚠ um gene não conta toda a história.

> **Fala:** "Resumindo: comparando as letras de um único gene, recuperamos os grupos que os biólogos conhecem, com um algoritmo que eu programei e validei. Os limites: um gene só, mitocondrial, herdado só da mãe, e um método simples que assume que todo mundo evolui na mesma velocidade. Estudos de verdade usam muitos genes. Mas a ideia central é essa: **a história da vida está escrita no DNA, e dá para ler com um pouco de código.** Obrigado! Perguntas?"

---

## Dicas e possíveis perguntas

- **Se o tempo apertar:** corte o slide 6 (validação) para uma frase e o slide 10 para 1 minuto.
- **"Por que gap = −2 e não −1?"** São parâmetros ajustáveis; usei match +1, mismatch −1, gap −2, que é uma escolha comum e penaliza buracos mais que trocas. Resultados mudam pouco no ranking de parentesco.
- **"Por que o mitocondrial?"** Existe em muitas cópias por célula (fácil de sequenciar), está presente em todos os vertebrados e muda em ritmo bom para comparar espécies.
- **"A baleia é parente do hipopótamo mesmo?"** Sim, os hipopótamos são os parentes vivos mais próximos das baleias (estudos com muitos genes confirmam). Nosso resultado empatado mostra o limite de usar um só gene.
- **"O que é UPGMA?"** Agrupar sempre os dois grupos mais próximos, usando a média das distâncias entre os membros. Simples, mas assume velocidade evolutiva igual.
- **Antes de apresentar:** abra `projeto_final.ipynb` no Colab, rode "Executar tudo" uma vez (leva ~1 min) e deixe a aba pronta como plano B para o caso de os slides falharem.
