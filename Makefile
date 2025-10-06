# ==== Compiler & Flags ====
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Wpedantic -IprofileLib -MMD -MP

# ==== Directories ====
SRC_DIR  := .
LIB_DIR  := profileLib
BUILD_DIR:= build


# ==== Sources & Targets ====
APP_NAME := mean_profile
TARGET   := $(APP_NAME)

SRCS := $(SRC_DIR)/main.cpp \
        $(LIB_DIR)/profile.cpp

# Turn e.g. profileLib/profile.cpp -> build/profileLib/profile.o
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# ==== Default rule ====
.PHONY: all
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@

# Compile .cpp -> .o into mirrored build/ tree
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Convenience targets ====
.PHONY: debug
debug: CXXFLAGS := -std=c++17 -g -O0 -Wall -Wextra -Wpedantic -IprofileLib -MMD -MP
debug: clean all

.PHONY: run
run: $(TARGET)
	./$(TARGET)

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR)

# Include auto-generated dep files
-include $(DEPS)
