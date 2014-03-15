CXX = g++
RM = rm -f
CPPFLAGS = -std=c++0x -g -MMD
LINK.o = g++

build_dir = build
test_dir = $(build_dir)/tests

target = filter
target_src = filter_main.cpp
target_obj = $(build_dir)/$(target_src:.cpp=.o)

srcs = $(wildcard *.cpp)
obj_srcs = $(filter-out $(target_src), $(srcs))
objs = $(addprefix $(build_dir)/, $(obj_srcs:.cpp=.o))
depends = $(srcs:.cpp=.d)

test_srcs = $(wildcard tests/*.cpp)
tests = $(addprefix $(build_dir)/, $(test_srcs:.cpp=))
test_depends = $(test_srcs:.cpp=.d)

vpath %.o $(build_dir)

.PHONY: all
all: $(build_dir) $(target)

.PHONY: clean
clean:
	$(RM) $(target) $(build_dir)/*.o $(build_dir)/*.d
	$(RM) $(test_dir)/*

$(target): $(target_obj) $(objs)
	$(CXX) $(CPPFLAGS) -o $@ $^ -lboost_system -lboost_thread-mt -lboost_program_options

$(build_dir):
	mkdir -p $@

.PRECIOUS: $(build_dir)/%.o
$(build_dir)/%.o: %.cpp
	$(CXX) $(CPPFLAGS) -o $@ -c $<

.PHONY: test
test: $(test_dir) $(tests)

$(test_dir):
	mkdir -p $@

$(test_dir)/%_test: $(test_dir)/%_test.o $(objs)
	$(CXX) $(CPPFLAGS) -o $@ $^
	$@

.PRECIOUS: $(test_dir)/%_test.o
$(test_dir)/%_test.o: tests/%_test.cpp %.hpp
	$(CXX) $(CPPFLAGS) -o $@ -c $<

-include $(depends)
-include $(test_depends)
