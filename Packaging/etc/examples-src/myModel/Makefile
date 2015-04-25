EXE=analysis

CXX=/usr/local/bin/mpicxx

CXXFLAGS=-Wall -O3 -DNDEBUG
CXXFLAGS+=`/Users/ayanpaul/Codes/Package/SusyFit-0.2/build/SusyFit/bin/susyfit-config --cflags` 
CXXFLAGS+=-I/Users/ayanpaul/Codes/Package/SusyFit-0.2/build/BAT/include -D_MPI -I/opt/local/include
CXXFLAGS+=`/usr/local/bin/root-config --cflags` 
CXXFLAGS+=`/usr/local/bin/gsl-config --cflags` 
CXXFLAGS+=-I ./src

LIBS=`/Users/ayanpaul/Codes/Package/SusyFit-0.2/build/SusyFit/bin/susyfit-config --libs`
LIBS+=-L/Users/ayanpaul/Codes/Package/SusyFit-0.2/build/BAT/lib -lBATmodels -lBAT  
LIBS+=`/usr/local/bin/root-config --libs` 
LIBS+=`/usr/local/bin/root-config --ldflags` 
LIBS+=`/usr/local/bin/gsl-config --libs` 
LIBS+=-Wl,-rpath,/Users/ayanpaul/Codes/Package/SusyFit-0.2/build/BAT/lib -Wl,-rpath,/usr/local/lib -Wl,-rpath,/usr/local/lib/root -Wl,-rpath,/opt/local/lib

all: $(EXE)

$(EXE): MCMC.o src/myModel.o src/gg4l.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm src/*.o *.o

distclean:
	rm $(EXE) src/*.o *.o
