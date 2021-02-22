# a build command to build randmst executable from randmst.cpp
all: myprogram.cpp
		g++ –o randmst randmst.cpp
clean:
		$(RM) randmst