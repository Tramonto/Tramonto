#
# a simple make file that compiles a set of java files
#

JAVAC=javac
JARC=jar
JAVA_FLAGS = -deprecation
IDEA_CLASSES=$(HOME)/src/Idea/Java/classes
CLASS_PATH = .:$(IDEA_CLASSES):$(IDEA_CLASSES)/idea.jar:$(IDEA_CLASSES)/xerces.jar

SOURCE = \
    Tramonto/TramontoDFTHandler.java \
    Tramonto/TramontoXMLParser.java \
    Tramonto/TramontoFluidParametersHandler.java \
    Tramonto/TramontoMixingParametersHandler.java \
    Tramonto/TramontoSurfaceParametersHandler.java \
    Tramonto/TramontoStateParametersHandler.java \
    Tramonto/TramontoComputationalParametersHandler.java 


OBJECTS = ${SOURCE:.java=.class}

.SUFFIXES: .java .class

.java.class :
	$(JAVAC) -classpath $(CLASS_PATH) $(JAVA_FLAGS) $<

all : $(OBJECTS)

jar : $(OBJECTS)
	jar cf Tramonto.jar $(OBJECTS)

clean :
	rm $(OBJECTS) *~

