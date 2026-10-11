/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * VERSION: v20.08.25
 * DESCRIPTION: Simple linear vector (generic type elements)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdbool.h>

#include "vector.h"

/*
 * Constants
 */
#define VECTOR_GROW_CAPACITY(cap) ((cap) ? ((cap) + ((cap) >> 1)) : 8)

/*
 * Error Checking
 */
static inline void vector_check_mul_u64(
    const uint64_t a,
    const uint64_t b) {
  if (VECTOR_UNLIKELY(b != 0 && a > UINT64_MAX/b)) {
    fprintf(stderr,"Vector allocation size overflow\n");
    exit(1);
  }
}

/*
 * Setup
 */
vector_t* vector_new_(
    const uint64_t num_initial_elements,
    const uint64_t element_size) {
  vector_t* const vector_buffer = (vector_t*) malloc(sizeof(vector_t));
  if (VECTOR_UNLIKELY(vector_buffer == NULL)) {
    fprintf(stderr,"Could not create new vector descriptor\n");
    exit(1);
  }
  vector_buffer->element_size = element_size;
  vector_buffer->elements_allocated = num_initial_elements;
  if (num_initial_elements == 0) {
    vector_buffer->memory = NULL;
  } else {
    vector_check_mul_u64(num_initial_elements,element_size);
    const uint64_t bytes = num_initial_elements*element_size;
    vector_buffer->memory = malloc(bytes);
    if (VECTOR_UNLIKELY(vector_buffer->memory == NULL)) {
      fprintf(stderr,"Could not create new vector (%" PRIu64 " bytes requested)",
          bytes);
      exit(1);
    }
  }
  vector_buffer->used = 0;
  return vector_buffer;
}
void vector_delete(
    vector_t* const vector) {
  free(vector->memory);
  free(vector);
}
void vector_cast(
    vector_t* const vector,
    const uint64_t element_size) {
  /*
   * Reuse the vector backing allocation for another element size.
   * Existing contents are discarded.
   */
  vector_check_mul_u64(vector->elements_allocated,vector->element_size);
  vector->elements_allocated = (vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size = element_size;
  vector->used = 0;
}
void vector_reserve(
    vector_t* const vector,
    const uint64_t num_elements,
    const bool zero_mem) {
  if (VECTOR_UNLIKELY(vector->elements_allocated < num_elements)) {
    uint64_t new_capacity = VECTOR_GROW_CAPACITY(vector->elements_allocated);
    if (new_capacity < num_elements) {
      new_capacity = num_elements;
    }
    vector_check_mul_u64(new_capacity,vector->element_size);
    const uint64_t bytes = new_capacity*vector->element_size;
    void* const new_memory = realloc(vector->memory,bytes);
    if (VECTOR_UNLIKELY(new_memory == NULL)) {
      fprintf(stderr,"Could not reserve vector (%" PRIu64 " bytes requested)",
          bytes);
      exit(1);
    }
    vector->memory = new_memory;
    vector->elements_allocated = new_capacity;
  }
  if (zero_mem && vector->elements_allocated > vector->used) {
    /*
     * zero_mem clears all currently unused capacity, not only newly allocated
     * bytes.
     */
    vector_check_mul_u64(vector->used,vector->element_size);
    vector_check_mul_u64(vector->elements_allocated-vector->used,vector->element_size);
    uint8_t* const base = (uint8_t*) vector->memory;
    memset(base+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
}
/*
 * Accessors
 */
#ifdef VECTOR_DEBUG
static void vector_check_element_size(
    vector_t* const vector,
    const uint64_t element_size) {
  if (VECTOR_UNLIKELY(element_size != vector->element_size)) {
    fprintf(stderr,
        "Vector element-size mismatch: requested %" PRIu64
        ", vector has %" PRIu64 "\n",
        element_size,vector->element_size);
    exit(1);
  }
}
void* vector_get_mem_element(
    vector_t* const vector,
    const uint64_t position,
    const uint64_t element_size) {
  vector_check_element_size(vector,element_size);
  if (position >= (vector)->used) {
    fprintf(stderr,"Vector position out-of-range [0,%"PRIu64")",(vector)->used);
    exit(1);
  }
  uint8_t* const base = (uint8_t*) vector->memory;
  return base + (position*element_size);
}
void* vector_get_mem_last_element(
    vector_t* const vector,
    const uint64_t element_size) {
  vector_check_element_size(vector,element_size);
  if (VECTOR_UNLIKELY(vector->used == 0)) {
    fprintf(stderr,"Vector last element requested from empty vector\n");
    exit(1);
  }
  uint8_t* const base = (uint8_t*) vector->memory;
  return base + ((vector->used-1)*element_size);
}
#endif
/*
 * Miscellaneous
 */
void vector_copy(
    vector_t* const vector_to,
    vector_t* const vector_from) {
  if (vector_to == vector_from) {
    return;
  }
  // Prepare
  vector_cast(vector_to,vector_from->element_size);
  vector_reserve(vector_to,vector_from->used,false);
  // Copy
  vector_set_used(vector_to,vector_from->used);
  vector_check_mul_u64(vector_from->used,vector_from->element_size);
  if (vector_from->used > 0) {
    memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
  }
}
vector_t* vector_dup(
    vector_t* const vector_src) {
  vector_t* const vector_cpy = vector_new_(vector_src->used,vector_src->element_size);
  // Copy
  vector_set_used(vector_cpy,vector_src->used);
  vector_check_mul_u64(vector_src->used,vector_src->element_size);
  if (vector_src->used > 0) {
    memcpy(vector_cpy->memory,vector_src->memory,vector_src->used*vector_src->element_size);
  }
  return vector_cpy;
}
