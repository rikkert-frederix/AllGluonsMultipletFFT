FC = gfortran
FFLAGS = -O3 -std=f2018 -Wall -Wextra
LDFLAGS =

BUILD_DIR = build
AMPLIGLUON_MULTIPLET_BUILD_DIR = build
AMPLIGLUON_TRACE_BUILD_DIR = build
WIGNER_TABLE ?= ../Wigner6j/data/su3_adjoint_swap_prefix_6.tbl

OBJECTS = \
	$(BUILD_DIR)/find_zero.o \
	$(BUILD_DIR)/rambo_on_diet.o

AMPLIGLUON_MULTIPLET_OBJECTS = \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/ampligluon_multiplet_kinds.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/wigner_table.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/multiplet_paths.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/recoupling_plan.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/multiplet_radiation.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/qcd_kinematics.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/multiplet_mhv.o \
	../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR)/ampligluon_multiplet.o

AMPLIGLUON_TRACE_LIBRARY = \
	../AmpliGluonTrace/$(AMPLIGLUON_TRACE_BUILD_DIR)/libampligluon_trace.a

.PHONY: all checker clean comparison_dependencies debug test

all: $(BUILD_DIR)/generate_ampligluon_events

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/find_zero.o: find_zero.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/rambo_on_diet.o: rambo_on_diet.f90 \
		$(BUILD_DIR)/find_zero.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/generate_ampligluon_events.o: generate_ampligluon_events.f90 \
		$(BUILD_DIR)/rambo_on_diet.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/generate_ampligluon_events: $(OBJECTS) \
		$(BUILD_DIR)/generate_ampligluon_events.o
	$(FC) $^ $(LDFLAGS) -o $@

comparison_dependencies:
	$(MAKE) -C ../AmpliGluonMultiplet \
		BUILD_DIR=$(AMPLIGLUON_MULTIPLET_BUILD_DIR) all
	$(MAKE) -C ../AmpliGluonTrace BUILD_DIR=$(AMPLIGLUON_TRACE_BUILD_DIR) all

$(BUILD_DIR)/event_input.o: event_input.f90 | $(BUILD_DIR) \
		comparison_dependencies
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) \
		-I../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR) \
		-c $< -o $@

$(BUILD_DIR)/compare_ampligluon_multiplet_trace.o: \
		compare_ampligluon_multiplet_trace.f90 \
		$(BUILD_DIR)/event_input.o | $(BUILD_DIR) \
		comparison_dependencies
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) \
		-I../AmpliGluonMultiplet/$(AMPLIGLUON_MULTIPLET_BUILD_DIR) \
		-I../AmpliGluonTrace/$(AMPLIGLUON_TRACE_BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/compare_ampligluon_multiplet_trace: \
		$(BUILD_DIR)/compare_ampligluon_multiplet_trace.o \
		$(BUILD_DIR)/event_input.o | \
		comparison_dependencies
	$(FC) $(AMPLIGLUON_MULTIPLET_OBJECTS) $^ $(AMPLIGLUON_TRACE_LIBRARY) \
		$(LDFLAGS) -o $@

checker: $(BUILD_DIR)/compare_ampligluon_multiplet_trace

test: $(BUILD_DIR)/generate_ampligluon_events \
		$(BUILD_DIR)/compare_ampligluon_multiplet_trace
	$(BUILD_DIR)/generate_ampligluon_events 2 2 \
		$(BUILD_DIR)/test_seed_a 1000 1729
	$(BUILD_DIR)/generate_ampligluon_events 2 2 \
		$(BUILD_DIR)/test_seed_b 1000 1729
	cmp $(BUILD_DIR)/test_seed_a_000001.event \
		$(BUILD_DIR)/test_seed_b_000001.event
	cmp $(BUILD_DIR)/test_seed_a_000002.event \
		$(BUILD_DIR)/test_seed_b_000002.event
	$(BUILD_DIR)/compare_ampligluon_multiplet_trace $(WIGNER_TABLE) \
		$(BUILD_DIR)/test_seed_a_000001.event
	$(BUILD_DIR)/compare_ampligluon_multiplet_trace $(WIGNER_TABLE) \
		$(BUILD_DIR)/test_seed_a_000002.event
	@echo 'RamboOnDiet event-generation regression: PASS'

debug:
	$(MAKE) BUILD_DIR=build-debug \
		AMPLIGLUON_MULTIPLET_BUILD_DIR=build-debug \
		AMPLIGLUON_TRACE_BUILD_DIR=build-debug \
		FFLAGS='-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all' test

clean:
	rm -rf build build-debug
