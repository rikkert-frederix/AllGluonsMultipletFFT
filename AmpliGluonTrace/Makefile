FC = gfortran
AR = ar
FFLAGS = -O3 -std=f2018 -Wall -Wextra
LDFLAGS =

BUILD_DIR = build

OBJECTS = \
	$(BUILD_DIR)/ampligluon_common.o \
	$(BUILD_DIR)/symmetric_group_fft.o \
	$(BUILD_DIR)/trace_colour_kernel.o \
	$(BUILD_DIR)/gluon_kinematics.o \
	$(BUILD_DIR)/trace_colour_matrix.o \
	$(BUILD_DIR)/trace_current_dag.o \
	$(BUILD_DIR)/trace_order_recursion.o \
	$(BUILD_DIR)/trace_mhv.o \
	$(BUILD_DIR)/ampligluon_trace.o

.PHONY: all clean debug test

all: $(BUILD_DIR)/libampligluon_trace.a

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/ampligluon_common.o: src/ampligluon_common.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gluon_kinematics.o: src/gluon_kinematics.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/symmetric_group_fft.o: src/symmetric_group_fft.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_colour_kernel.o: src/trace_colour_kernel.f90 \
		$(BUILD_DIR)/ampligluon_common.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_colour_matrix.o: src/trace_colour_matrix.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/symmetric_group_fft.o \
		$(BUILD_DIR)/trace_colour_kernel.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_current_dag.o: src/trace_current_dag.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_order_recursion.o: src/trace_order_recursion.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/trace_mhv.o: src/trace_mhv.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/trace_order_recursion.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/ampligluon_trace.o: src/ampligluon_trace.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o \
		$(BUILD_DIR)/trace_colour_matrix.o \
		$(BUILD_DIR)/trace_current_dag.o \
		$(BUILD_DIR)/trace_order_recursion.o \
		$(BUILD_DIR)/trace_mhv.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/libampligluon_trace.a: $(OBJECTS)
	$(AR) rcs $@ $^

$(BUILD_DIR)/test_ampligluon_trace.o: tests/test_ampligluon_trace.f90 \
		$(BUILD_DIR)/ampligluon_trace.o \
		$(BUILD_DIR)/trace_colour_kernel.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_ampligluon_trace: \
		$(BUILD_DIR)/test_ampligluon_trace.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_symmetric_group_fft.o: \
		tests/test_symmetric_group_fft.f90 \
		$(BUILD_DIR)/symmetric_group_fft.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_symmetric_group_fft: \
		$(BUILD_DIR)/test_symmetric_group_fft.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_trace_colour_kernel.o: \
		tests/test_trace_colour_kernel.f90 \
		$(BUILD_DIR)/trace_colour_kernel.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_trace_colour_kernel: \
		$(BUILD_DIR)/test_trace_colour_kernel.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_trace_colour_matrix.o: \
		tests/test_trace_colour_matrix.f90 \
		$(BUILD_DIR)/trace_colour_matrix.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_trace_colour_matrix: \
		$(BUILD_DIR)/test_trace_colour_matrix.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_trace_current_dag.o: \
		tests/test_trace_current_dag.f90 \
		$(BUILD_DIR)/trace_current_dag.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_trace_current_dag: \
		$(BUILD_DIR)/test_trace_current_dag.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/legacy_trace_order.o: tests/legacy_trace_order.f90 \
		$(BUILD_DIR)/ampligluon_common.o \
		$(BUILD_DIR)/gluon_kinematics.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_trace_current_equivalence.o: \
		tests/test_trace_current_equivalence.f90 \
		$(BUILD_DIR)/legacy_trace_order.o \
		$(BUILD_DIR)/trace_current_dag.o \
		$(BUILD_DIR)/trace_order_recursion.o \
		$(BUILD_DIR)/trace_mhv.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_trace_current_equivalence: \
		$(BUILD_DIR)/test_trace_current_equivalence.o \
		$(BUILD_DIR)/legacy_trace_order.o \
		$(BUILD_DIR)/libampligluon_trace.a
	$(FC) $^ $(LDFLAGS) -o $@

test: $(BUILD_DIR)/test_ampligluon_trace \
		$(BUILD_DIR)/test_symmetric_group_fft \
		$(BUILD_DIR)/test_trace_colour_kernel \
		$(BUILD_DIR)/test_trace_colour_matrix \
		$(BUILD_DIR)/test_trace_current_dag \
		$(BUILD_DIR)/test_trace_current_equivalence
	$(BUILD_DIR)/test_ampligluon_trace
	$(BUILD_DIR)/test_symmetric_group_fft
	$(BUILD_DIR)/test_trace_colour_kernel
	$(BUILD_DIR)/test_trace_colour_matrix
	$(BUILD_DIR)/test_trace_current_dag
	$(BUILD_DIR)/test_trace_current_equivalence

debug:
	$(MAKE) BUILD_DIR=build-debug \
		FFLAGS='-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all' test

clean:
	rm -rf build build-debug
