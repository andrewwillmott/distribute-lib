TARGET = distribute
CXXFLAGS = --std=c++11

all: $(TARGET)

$(TARGET): Distribute.cpp Distribute.hpp Generate.hpp DistributeTest.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) Distribute.cpp Generate.cpp DistributeTest.cpp

clean:
	$(RM) $(TARGET)

test: distribute
	./distribute > dist-test.txt
	diff dist-test.txt dist-test-ref.txt