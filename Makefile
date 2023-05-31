PKG_CONFIG ?= pkg-config
DEPS = gmock gtest_main
LIBS += $(shell $(PKG_CONFIG) --libs $(DEPS))
CXXFLAGS += $(shell $(PKG_CONFIG) --cflags $(DEPS))
CXXFLAGS += -Wall -Wextra -pedantic -std=c++20
DEST = gf256_test
SOURCES = gf256_test.cc
HEADERS = gf256.h


all: $(DEST)
	./$(DEST)

$(DEST): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(SOURCES) -o $@

clean:
	rm -f $(DEST)

.PHONY: all clean
