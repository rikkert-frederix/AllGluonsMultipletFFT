FC = gfortran
AR = ar
FFLAGS = -O3 -std=f2018 -Wall -Wextra
LDFLAGS =

BUILD_DIR = build
AMPLIGLUON_TRACE_SOURCE_DIR = ../AmpliGluonTrace/src

OBJECTS = \
	$(BUILD_DIR)/ampligluon_common.o \
	$(BUILD_DIR)/symmetric_group_fft.o \
	$(BUILD_DIR)/trace_colour_kernel.o \
	$(BUILD_DIR)/gluon_kinematics.o \
	$(BUILD_DIR)/trace_current_dag.o \
	$(BUILD_DIR)/adjoint_current_dag.o \
	$(BUILD_DIR)/adjoint_colour_matrix.o \
	$(BUILD_DIR)/ampligluon_adjoint.o

.PHONY: all clean debug test

all: $(BUILD_DIR)/libampligluon_adjoint.a

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/ampligluon_common.o: \
		$(AMPLIGLUON_TRACE_SOURCE_DIR)/ampligluon_common.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gluon_kinematics.o: \
		$(AMPLIGLUON_TRACE_SOURCE_DIR)/gluon_kinematics.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/symmetric_group_fft.o: \
		$(AMPLIGLUON_TRACE_SOURCE_DIR)/symmetric_group_fft.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_colour_kernel.o: \
		$(AMPLIGLUON_TRACE_SOURCE_DIR)/trace_colour_kernel.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/adjoint_colour_matrix.o: src/adjoint_colour_matrix.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/symmetric_group_fft.o \
		$(BUILD_DIR)/trace_colour_kernel.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_current_dag.o: \
		$(AMPLIGLUON_TRACE_SOURCE_DIR)/trace_current_dag.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/adjoint_current_dag.o: src/adjoint_current_dag.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/trace_current_dag.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/ampligluon_adjoint.o: src/ampligluon_adjoint.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o \
		$(BUILD_DIR)/adjoint_current_dag.o \
		$(BUILD_DIR)/adjoint_colour_matrix.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/libampligluon_adjoint.a: $(OBJECTS)
	$(AR) rcs $@ $^

$(BUILD_DIR)/test_ampligluon_adjoint.o: tests/test_ampligluon_adjoint.f90 \
		$(BUILD_DIR)/ampligluon_adjoint.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_ampligluon_adjoint: \
		$(BUILD_DIR)/test_ampligluon_adjoint.o \
		$(BUILD_DIR)/libampligluon_adjoint.a
	$(FC) $^ $(LDFLAGS) -o $@

test: $(BUILD_DIR)/test_ampligluon_adjoint
	$(BUILD_DIR)/test_ampligluon_adjoint

debug:
	$(MAKE) BUILD_DIR=build-debug \
		FFLAGS='-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all' test

clean:
	rm -rf build build-debug
