/*------------------------------------------------------------------------*/
/*! \file aig.cpp
    \brief contains functions to parse and manipulate the input AIG

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "aig.h"
#include "gate.h"  // signed_mult
/*------------------------------------------------------------------------*/
// Global Variables
unsigned M, NN;
unsigned a0, al, ainc;
unsigned b0, bl, binc;

aiger * rewritten;
aiger * miter;
/*------------------------------------------------------------------------*/
// Local Variables

aiger * model;   // /< aiger* object, used for storing the given AIG graph
/*------------------------------------------------------------------------*/
void init_aig_substitution() {
  rewritten = aiger_init();
  miter = aiger_init();
}
/*------------------------------------------------------------------------*/
void reset_aig_substitution() {
  aiger_reset(rewritten);
  aiger_reset(miter);
}
/*------------------------------------------------------------------------*/
void init_aig_parsing() {
  model = aiger_init();
}
/*------------------------------------------------------------------------*/
void reset_aig_parsing() {
  aiger_reset(model);
}
/*------------------------------------------------------------------------*/
const char * aiger_open_and_read_to_model(const char * input_name) {
  return aiger_open_and_read_from_file(model, input_name);
}
/*------------------------------------------------------------------------*/

bool is_model_input(unsigned val) {
  return aiger_is_input(model, val);
}
/*------------------------------------------------------------------------*/
aiger_and * is_model_and(unsigned val) {
  return aiger_is_and(model, val);
}
/*------------------------------------------------------------------------*/
unsigned get_model_num_latches() {return model->num_latches;}
unsigned get_model_num_inputs() {return model->num_inputs;}
unsigned get_model_num_ands()   {return model->num_ands;}
unsigned get_model_num_outputs() {return model->num_outputs;}
unsigned get_model_maxvar()     {return model->maxvar;}
/*------------------------------------------------------------------------*/
unsigned get_model_inputs_lit(unsigned i) {
  assert(i < NN);
  return model->inputs[i].lit;
}
/*------------------------------------------------------------------------*/
const char* get_model_inputs_name(unsigned i) {
  assert(i < NN);
  return model->inputs[i].name;
}
/*------------------------------------------------------------------------*/
// Random-simulation bug check of the input AIG against the multiplier
// specification (a*b mod 2^NN, two's complement if -signed). Used before the
// Booth guard re-runs the algebra with another substitution order: if a
// differing input vector is found, the caller reports INCORRECT MULTIPLIER
// with that vector as counter example. 64 vectors per round (bit-parallel),
// splitmix64 inputs; only for NN/2 <= 64.
// Returns 1 and fills cex (one character per input, in input order).
const char * input_aig_name = 0;
static unsigned long long splitmix64_next(unsigned long long & x) {
  unsigned long long z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
int sim_check_multiplier(unsigned rounds, unsigned long long seed,
                         std::vector<char> & cex) {
  if (!model || NN/2 > 64 || model->num_outputs < NN) return 0;
  const unsigned n = NN/2;
  std::vector<unsigned long long> val(model->maxvar + 1, 0);
  auto lv = [&](unsigned l) -> unsigned long long {
    unsigned long long v = val[l >> 1];
    return (l & 1) ? ~v : v;
  };
  for (unsigned r = 0; r < rounds; r++) {
    for (unsigned j = 0; j < model->num_inputs; j++)
      val[model->inputs[j].lit >> 1] = splitmix64_next(seed);
    for (unsigned k = 0; k < model->num_ands; k++) {
      aiger_and & g = model->ands[k];
      val[g.lhs >> 1] = lv(g.rhs0) & lv(g.rhs1);
    }
    for (unsigned lane = 0; lane < 64; lane++) {
      unsigned __int128 a = 0, b = 0;
      for (unsigned i = 0; i < n; i++) {
        if ((val[model->inputs[a0 + i*ainc].lit >> 1] >> lane) & 1)
          a |= (unsigned __int128) 1 << i;
        if ((val[model->inputs[b0 + i*binc].lit >> 1] >> lane) & 1)
          b |= (unsigned __int128) 1 << i;
      }
      if (signed_mult && n < 128) {
        if ((a >> (n-1)) & 1) a |= ~(unsigned __int128) 0 << n;
        if ((b >> (n-1)) & 1) b |= ~(unsigned __int128) 0 << n;
      }
      unsigned __int128 pr = a * b;
      for (unsigned k = 0; k < NN; k++) {
        unsigned o = (lv(model->outputs[k].lit) >> lane) & 1;
        if (o != (unsigned) ((pr >> k) & 1)) {
          cex.assign(model->num_inputs, '0');
          for (unsigned j = 0; j < model->num_inputs; j++)
            if ((val[model->inputs[j].lit >> 1] >> lane) & 1) cex[j] = '1';
          return 1;
        }
      }
    }
  }
  return 0;
}
/*------------------------------------------------------------------------*/
unsigned slit(unsigned i) {
  if (!model) return -1;
  assert(i < NN);
  return model->outputs[i].lit;
}
/*------------------------------------------------------------------------*/
int write_model(FILE *file) {
  return aiger_write_to_file(model, aiger_binary_mode, file);
}
/*------------------------------------------------------------------------*/
