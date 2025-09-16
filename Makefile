CXXFLAGS = --std=c++11

DIST_INCLUDES=Generate.hpp Distribute.hpp
DIST_SOURCES=Generate.cpp Distribute.cpp

distribute: $(DIST_INCLUDES) $(DIST_SOURCES) DistributeTest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $(DIST_SOURCES) DistributeTest.cpp

clean:
	$(RM) distribute

test: distribute
	./distribute > dist-test.txt
	diff dist-test.txt dist-test-ref.txt