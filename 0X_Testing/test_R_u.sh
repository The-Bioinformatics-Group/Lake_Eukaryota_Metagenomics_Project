dockR = docker run -it -w $PWD --rm -v /home/$USER:/home/$USER dada2:1.6.0-1-gbadb111 R

dockR --no-save test_R.R > test_R.out


