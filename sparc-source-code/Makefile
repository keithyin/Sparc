# ===== toolchain =====
CXX := g++
AR  := ar

# ===== flags =====
CXXFLAGS := -O3 -std=c++14 -fPIC

# ===== target =====
TARGET := libsparc.a

# ===== sources =====
SRCS := $(filter-out main.cpp,$(wildcard *.cpp))
OBJS := $(SRCS:.cpp=.o)

# ===== rules =====
all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) rcs $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
