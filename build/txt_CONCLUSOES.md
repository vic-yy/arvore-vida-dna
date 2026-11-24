## 6. Conclusões e limitações

1. Comparar o DNA funciona como "régua de parentesco": quanto mais parecidos os genes, mais próximos os animais. Nossa árvore reencontrou os grupos conhecidos: primatas, felinos, cães, e os mamíferos separados da galinha e do peixe.
2. O algoritmo de Needleman-Wunsch, implementado do zero, deu **exatamente a mesma pontuação** que o Biopython em todos os pares testados, e o NumPy o tornou muito mais rápido, mantendo o resultado idêntico.
3. O tempo cresce como n·m: dobrar o tamanho das sequências quadruplica o tempo, o que limita o uso em dados gigantes.
4. **Limitações:** um único gene (CYTB) não conta toda a história evolutiva; ele é do DNA mitocondrial, herdado só da mãe, e muda em ritmos diferentes em cada linhagem (o camundongo e o par hipopótamo-baleia ficaram fora do esperado).
5. Além disso, usamos só a identidade simples (sem corrigir mutações repetidas) e o UPGMA, que assume que todas as linhagens evoluem na mesma velocidade. Estudos sérios usam muitos genes e métodos mais sofisticados.
6. Sobre os dados: para humano, gorila, camundongo e peixe-zebra o CYTB foi extraído de genomas mitocondriais de referência; os demais vêm de registros do gene, com 1.140 a 1.144 bases em todas as espécies.