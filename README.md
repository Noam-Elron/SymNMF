#TODO:
1. Analysis of KMeans vs Symnmf.

Useful Commands:
1. Valgrind Python Memory Leak Checking: valgrind python3 --suppressions=/usr/lib/valgrind/python3.supp symnmf.py 2 symnmf ./tests/simple_test2.txt
2. Valgrind C Direct Interface Memory Leak Checking: valgrind --leak-check=full ./test ddg ./tests/simple_test2.txt 
3. Python C-Extension Compilation: python3 setup.py build_ext --inplace
4. C Direct Interface Compilation: gcc -ansi -Wall -Wextra -Werror -pedantic-errors utils.c sym.c norm.c diagonal.c symnmf.c -o test -lm
