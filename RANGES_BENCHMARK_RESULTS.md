# Ranges Performance Benchmark Results

## Overview

This document contains the performance analysis of the current `Ranges` implementation in fibertools-rs, including benchmarks and optimization recommendations.

## Benchmark Results

### 1. Ranges Creation Performance

**Scale**: Creating Ranges with different numbers of features (1-1000)

| Features | Time (µs) | Notes |
|----------|-----------|-------|
| 1        | 0.862     | Excellent for single features |
| 5        | 1.611     | Good scalability |
| 10       | 2.481     | Linear scaling continues |
| 50       | 7.380     | Still reasonable |
| 100      | 15.354    | Starting to show overhead |
| 500      | 88.441    | Significant cost increase |
| 1000     | 181.25    | **Bottleneck identified** |

**Key Insights**:
- Creation time scales roughly linearly with feature count
- Overhead becomes significant at 500+ features
- Main cost is likely from coordinate transformation and vector allocations

### 2. Getter Method Performance

**Scale**: Accessing data from 100-feature Ranges

| Method               | Time (ps/ns) | Performance Level |
|---------------------|--------------|------------------|
| `starts()`          | 632 ps       | ⚡ Excellent (direct reference) |
| `ends()`            | 662 ps       | ⚡ Excellent (direct reference) |
| `lengths()`         | 661 ps       | ⚡ Excellent (direct reference) |
| `reference_starts()` | 663 ps       | ⚡ Excellent (direct reference) |
| `get_starts()`      | 258 ns       | ⚠️ 400x slower (clone + flatten) |
| `get_molecular()`   | 66 ns        | ✅ Good (iterator-based) |
| `get_reference()`   | 71 ns        | ✅ Good (iterator-based) |

**Key Insights**:
- Direct reference getters are nearly free (~650 picoseconds)
- `get_starts()` is expensive due to cloning and flattening
- Iterator-based methods are reasonable compromise

### 3. Iteration Performance

**Scale**: Iterating over different numbers of features

| Features | Time (µs) | Per-feature (ns) |
|----------|-----------|------------------|
| 10       | 0.021     | 2.1             |
| 100      | 0.214     | 2.1             |
| 1000     | 2.160     | 2.2             |

**Key Insights**:
- Iteration scales perfectly linearly
- ~2.1ns per feature is excellent performance
- Iterator implementation is well-optimized

### 4. Filtering Performance

**Scale**: Filtering 500-feature Ranges

| Operation                     | Time (µs) | Notes |
|-------------------------------|-----------|-------|
| `filter_by_qual()`           | 89.9      | Rebuilds all vectors |
| `filter_starts_at_read_ends()` | 89.9      | Similar rebuild cost |

**Key Insights**:
- Filtering is expensive due to vector rebuilding
- Both operations have similar performance profiles
- Opportunity for optimization with better data structures

### 5. Cloning Performance

**Scale**: Cloning Ranges with different feature counts

| Features | Time (µs) | Memory Impact |
|----------|-----------|---------------|
| 10       | 0.206     | Minimal |
| 100      | 0.254     | Low |
| 1000     | 1.420     | Moderate |

**Key Insights**:
- Cloning scales sub-linearly (good!)
- Even 1000 features clone in <1.5µs
- Memory allocation is well-optimized

### 6. String Serialization Performance

**Scale**: Converting 100-feature Ranges to strings

| Operation              | Time (µs) | Notes |
|------------------------|-----------|-------|
| `to_strings_reference` | 9.20      | Slightly faster |
| `to_strings_molecular` | 11.17     | Includes quality data |

**Key Insights**:
- String conversion is moderately expensive
- Quality data adds ~20% overhead
- Opportunity for streaming/buffered output

### 7. Merge Performance

**Scale**: Merging multiple Ranges

| Ranges Count | Time (µs) | Notes |
|--------------|-----------|-------|
| 2            | 5.01      | Base case |
| 5            | 20.31     | 4x increase |
| 10           | 46.17     | 2.3x increase |
| 20           | 137.92    | 3x increase |

**Key Insights**:
- Merge performance degrades super-linearly
- Sorting step is likely the bottleneck
- Could benefit from more efficient merge algorithms

### 8. Memory Pattern Analysis

**Scale**: Allocation patterns

| Pattern            | Time (ms) | Notes |
|--------------------|-----------|-------|
| Many small ranges  | 0.319     | 100 ranges × 5 features |
| Few large ranges   | 1.064     | 10 ranges × 500 features |

**Key Insights**:
- Small ranges are more memory-efficient
- Large ranges have higher overhead per feature
- Memory fragmentation may be an issue

## Performance Bottlenecks Identified

### 1. **Critical**: Large Ranges Creation (500+ features)
- **Impact**: 88-181µs for 500-1000 features
- **Cause**: Coordinate transformation overhead
- **Fix Priority**: High

### 2. **Major**: `get_starts()` Method
- **Impact**: 400x slower than direct access
- **Cause**: Clone + flatten operations
- **Fix Priority**: High

### 3. **Moderate**: Filtering Operations
- **Impact**: ~90µs for 500 features
- **Cause**: Full vector rebuilding
- **Fix Priority**: Medium

### 4. **Moderate**: Merge Operations
- **Impact**: Super-linear scaling
- **Cause**: Inefficient sorting/merging
- **Fix Priority**: Medium

### 5. **Minor**: String Serialization
- **Impact**: ~10µs for 100 features
- **Cause**: String allocations
- **Fix Priority**: Low

## Optimization Recommendations

### Immediate (High Impact, Low Effort)

1. **Replace `get_starts()` with iterator**
   ```rust
   // Instead of:
   pub fn get_starts(&self) -> Vec<i64> {
       self.starts.clone().into_iter().flatten().collect()
   }
   
   // Use:
   pub fn starts_iter(&self) -> impl Iterator<Item = i64> + '_ {
       self.starts.iter().filter_map(|&x| x)
   }
   ```

2. **Add `smallvec` for small collections**
   ```rust
   use smallvec::SmallVec;
   type RangeVec<T> = SmallVec<[T; 8]>;  // Stack allocation for ≤8 elements
   ```

### Medium-term (High Impact, Medium Effort)

3. **Struct-of-Arrays to Array-of-Structs**
   ```rust
   #[derive(Clone, Debug)]
   pub struct Ranges {
       features: Vec<Feature>,
       seq_len: i64,
       reverse: bool,
   }
   
   #[derive(Clone, Debug)]
   struct Feature {
       start: Option<i64>,      // Keep i64 for genomic coordinates
       length: Option<i64>,     // derive end from start+length
       qual: u8,
       reference: Option<ReferenceFeature>,
   }
   ```

4. **Lazy Reference Coordinate Calculation**
   ```rust
   // Only compute reference coordinates when needed
   reference_cache: OnceCell<Vec<Option<ReferenceFeature>>>,
   ```

5. **Optimize Filtering with Bit Vectors**
   ```rust
   // Use bit vectors for keep/drop decisions
   use bit_vec::BitVec;
   
   fn filter_with_bitvec(&mut self, keep: &BitVec) {
       // More efficient than rebuilding vectors
   }
   ```

### Long-term (Medium Impact, High Effort)

6. **Memory Pool for Frequent Allocations**
   ```rust
   // Pre-allocate common sizes
   thread_local! {
       static RANGE_POOL: RefCell<Vec<Vec<Feature>>> = RefCell::new(Vec::new());
   }
   ```

7. **SIMD Optimizations for Coordinate Transformations**
   ```rust
   // Use SIMD for bulk coordinate transformations
   #[cfg(target_arch = "x86_64")]
   use std::arch::x86_64::*;
   ```

## Expected Performance Gains

| Optimization | Expected Improvement | Effort |
|--------------|---------------------|--------|
| Iterator-based getters | 300-400x for `get_starts()` | Low |
| SmallVec | 20-50% for small ranges | Low |
| Remove redundant fields | 25-30% memory, 10-15% speed | Medium |
| Struct-of-arrays restructure | 15-25% better cache locality | Medium |
| Lazy reference calc | 10-30% for creation | Medium |
| Filtering optimization | 50-80% for filter ops | High |

## Conclusion

The current `Ranges` implementation performs well for small to medium datasets (≤100 features) but shows performance bottlenecks at larger scales. The most impactful optimizations would be:

1. **Immediate**: Replace cloning getters with iterators
2. **Short-term**: Add `smallvec` for memory efficiency  
3. **Medium-term**: Restructure to remove redundancy and improve cache locality

These changes could provide 2-5x performance improvements for common use cases while maintaining API compatibility.