# Executa um notebook: python run_nb.py entrada.ipynb saida.ipynb  (cwd = pasta do projeto)
import sys, nbformat
from nbconvert.preprocessors import ExecutePreprocessor
nb = nbformat.read(sys.argv[1], as_version=4)
ExecutePreprocessor(timeout=1800, kernel_name="python3").preprocess(nb, {"metadata": {"path": "."}})
nbformat.write(nb, sys.argv[2])
print("executado ->", sys.argv[2])
