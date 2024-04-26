all: compile

compile:
	mpicxx -std=c++20 life.cpp -o life

sync:
	rsync --exclude 'build' -rvu ~/prl/projekt2 xkrato61@merlin.fit.vutbr.cz:~/prl/