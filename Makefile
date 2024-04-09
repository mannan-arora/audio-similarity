CXX = g++
CXXFLAGS = -std=c++20 -g -Wall -MMD -Werror=vla
OBJECTS = main.o similarity.o
DEPENDS = ${OBJECTS:.o=.d}
EXEC = similarity

${EXEC} : ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o ${EXEC} -lsndfile -lfftw3

clean :
	rm -f ${DEPENDS} ${OBJECTS} ${EXEC}

-include ${DEPENDS} # reads the .d files and reruns dependencies
