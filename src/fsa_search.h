/*------------------------------------------------------------------------*/
/*! \file fsa_search.h
    \brief contains a search-based detection of the final-stage adder

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
  Copyright(C) 2026 Ben Jensen (final-stage-adder search)
  Released under the MIT License, see LICENSE.
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_FSA_SEARCH_H_
#define AMULET2_SRC_FSA_SEARCH_H_
#include <stdio.h>
#include "aig.h"
/**
    Replaces identify_final_stage_adder() + build_adder_miter():
    finds, column by column, the two rows (x_i, y_i) entering a final-stage
    adder such that the outputs s_0..s_{n-1} equal the ripple sum of the rows
    (checked by simulation of the cone between rows and outputs on random row
    values), then builds the adder miter (region vs ripple-carry adder over the
    same rows) in 'miter' and the rewritten AIG (rows' cones + ripple adder)
    in 'rewritten'.
    @return true if an adder with at least one two-row column was found
*/
bool fsa_search_substitute();
/** Writes the adder miter as DIMACS; the header counts the clauses written. */
bool fsa_search_miter_to_file(FILE * file);
/** Writes the rewritten AIG prepared by fsa_search_substitute(). */
bool fsa_search_write_rewritten(FILE * file);
#endif  // AMULET2_SRC_FSA_SEARCH_H_
