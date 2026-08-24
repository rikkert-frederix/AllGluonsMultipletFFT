FC = gfortran
FFLAGS = -O3 -std=f2018 -Wall -Wextra
LDFLAGS =

BUILD_DIR = build
WIGNER_TABLE ?= ../Wigner6j/data/su3_adjoint_swap_prefix_4.tbl
SHALLOW_WIGNER_TABLE ?= ../Wigner6j/data/su3_adjoint_swap_prefix_3.tbl
MINIMAL_WIGNER_TABLE ?= ../Wigner6j/data/su3_adjoint_swap_prefix_2.tbl

COLOUR_OBJECTS = \
	$(BUILD_DIR)/ampligluon_multiplet_kinds.o \
	$(BUILD_DIR)/wigner_table.o \
	$(BUILD_DIR)/multiplet_paths.o \
	$(BUILD_DIR)/recoupling_plan.o \
	$(BUILD_DIR)/multiplet_radiation.o

CORE_OBJECTS = \
	$(COLOUR_OBJECTS) \
	$(BUILD_DIR)/qcd_kinematics.o \
	$(BUILD_DIR)/multiplet_mhv.o \
	$(BUILD_DIR)/ampligluon_multiplet.o

.PHONY: all clean debug native test

all: $(BUILD_DIR)/test_recoupling $(BUILD_DIR)/test_multiplet_radiation \
		$(BUILD_DIR)/test_ampligluon_multiplet

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/ampligluon_multiplet_kinds.o: \
		src/ampligluon_multiplet_kinds.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/wigner_table.o: src/wigner_table.f90 \
		$(BUILD_DIR)/ampligluon_multiplet_kinds.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/multiplet_paths.o: src/multiplet_paths.f90 \
		$(BUILD_DIR)/ampligluon_multiplet_kinds.o \
		$(BUILD_DIR)/wigner_table.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/recoupling_plan.o: src/recoupling_plan.f90 \
		$(BUILD_DIR)/ampligluon_multiplet_kinds.o \
		$(BUILD_DIR)/wigner_table.o \
		$(BUILD_DIR)/multiplet_paths.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/multiplet_radiation.o: src/multiplet_radiation.f90 \
		$(BUILD_DIR)/ampligluon_multiplet_kinds.o \
		$(BUILD_DIR)/multiplet_paths.o \
		$(BUILD_DIR)/recoupling_plan.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/qcd_kinematics.o: src/qcd_kinematics.f90 \
		$(BUILD_DIR)/ampligluon_multiplet_kinds.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/multiplet_mhv.o: src/multiplet_mhv.f90 \
		$(COLOUR_OBJECTS) $(BUILD_DIR)/qcd_kinematics.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/ampligluon_multiplet.o: src/ampligluon_multiplet.f90 \
		$(COLOUR_OBJECTS) $(BUILD_DIR)/qcd_kinematics.o \
		$(BUILD_DIR)/multiplet_mhv.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_recoupling.o: tests/test_recoupling.f90 \
		$(COLOUR_OBJECTS) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_ampligluon_multiplet.o: \
		tests/test_ampligluon_multiplet.f90 \
		$(CORE_OBJECTS) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_multiplet_radiation.o: \
		tests/test_multiplet_radiation.f90 \
		$(COLOUR_OBJECTS) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_recoupling: $(COLOUR_OBJECTS) \
		$(BUILD_DIR)/test_recoupling.o
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_ampligluon_multiplet: $(CORE_OBJECTS) \
		$(BUILD_DIR)/test_ampligluon_multiplet.o
	$(FC) $^ $(LDFLAGS) -o $@

$(BUILD_DIR)/test_multiplet_radiation: $(COLOUR_OBJECTS) \
		$(BUILD_DIR)/test_multiplet_radiation.o
	$(FC) $^ $(LDFLAGS) -o $@

test: $(BUILD_DIR)/test_recoupling \
		$(BUILD_DIR)/test_multiplet_radiation \
		$(BUILD_DIR)/test_ampligluon_multiplet
	$(BUILD_DIR)/test_recoupling $(WIGNER_TABLE)
	$(BUILD_DIR)/test_multiplet_radiation $(WIGNER_TABLE)
	$(BUILD_DIR)/test_ampligluon_multiplet $(WIGNER_TABLE)
	$(BUILD_DIR)/test_ampligluon_multiplet $(MINIMAL_WIGNER_TABLE) 5
	$(BUILD_DIR)/test_ampligluon_multiplet $(SHALLOW_WIGNER_TABLE) 7
	@if $(BUILD_DIR)/test_ampligluon_multiplet $(SHALLOW_WIGNER_TABLE) 8; then \
		echo 'Wigner table-depth boundary regression: FAIL'; \
		exit 1; \
	else \
		echo 'Wigner table-depth boundary regression: PASS'; \
	fi
	@echo 'table-depth regression: PASS'

debug:
	$(MAKE) BUILD_DIR=build-debug \
		FFLAGS='-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all' all test

native:
	$(MAKE) BUILD_DIR=build-native \
		FFLAGS='-Ofast -march=native -std=f2018 -Wall -Wextra' all

clean:
	rm -rf build build-debug build-native
