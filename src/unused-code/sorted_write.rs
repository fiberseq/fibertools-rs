use std::cmp::Reverse;
use std::collections::BinaryHeap;

pub trait SortedWriterTrait {
    fn get_chrom(&self) -> &str;
}
pub struct SortedWriter<T> {
    pub heap: BinaryHeap<Reverse<T>>,
    pub first_of_next: Option<Reverse<T>>,
}

impl<T: Ord + Clone> SortedWriter<T> {
    pub fn new() -> Self {
        Self {
            heap: BinaryHeap::new(),
            first_of_next: None,
        }
    }

    pub fn add_next_and_get_available(
        &mut self,
        mut values: Vec<T>,
        new_next: T,
    ) -> Option<Vec<T>> {
        // should be very fast (linear?) if presorted
        values.sort();

        self.first_of_next = Some(Reverse(new_next));
        for v in values {
            self.heap.push(Reverse(v));
        }
        self.pop()
    }

    /// return values that are smaller than `start_of_next`
    fn pop(&mut self) -> Option<Vec<T>> {
        let first_of_next = self.first_of_next.as_ref()?;
        let mut rtn = vec![];
        while self.heap.is_empty() == false {
            let min_val = self.heap.peek()?;
            if min_val.0 < first_of_next.0 {
                let owned_min_val = self.heap.pop()?;
                rtn.push(owned_min_val.0);
            } else {
                break;
            }
        }
        Some(rtn)
    }

    pub fn flush(self) -> Option<Vec<T>> {
        let rtn: Vec<T> = self
            .heap
            .into_sorted_vec()
            .into_iter()
            .rev()
            .map(|v| v.0)
            .collect();
        Some(rtn)
    }
}

/// a set of test cases for the SortedWriter
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let mut sw = SortedWriter::new();
        let mut rtn = sw.add_next_and_get_available(vec![1, 2, 3], 0);
        assert_eq!(rtn, Some(vec![]));
        rtn = sw.add_next_and_get_available(vec![4, 5, 6], 3);
        assert_eq!(rtn, Some(vec![1, 2]));
        rtn = sw.add_next_and_get_available(vec![4, 7, 8, 9], 4);
        assert_eq!(rtn, Some(vec![3]));
        rtn = sw.add_next_and_get_available(vec![10, 11, 12], 10);
        assert_eq!(rtn, Some(vec![4, 4, 5, 6, 7, 8, 9]));
        rtn = sw.flush();
        assert_eq!(rtn, Some(vec![10, 11, 12]));
    }
}
