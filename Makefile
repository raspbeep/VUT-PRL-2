CV_DIR=cv2

make: build
	cp ./build/gol ./gol

sync:
	# rsync --exclude 'build' -rvu ~/ppp/${CV_DIR}/ xkrato61@barbora.it4i.cz:~/ppp/${CV_DIR}
	# rsync --exclude 'build' -rvu ~/ppp/${CV_DIR}/ xkrato61@merlin.fit.vutbr.cz:~/ppp/${CV_DIR}
	rsync --exclude 'build' -rvu ~/ppp/${CV_DIR}/ xkrato61@barbora.it4i.cz:~/ppp/${CV_DIR}

gen:
	cmake -Bbuild -S.

build:
	cmake --build build --config Release