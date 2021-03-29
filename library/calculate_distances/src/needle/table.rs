//! Two-dimensional table utilities.
use core::fmt::{self, Debug};
use core::ops::Index;
use core::ops::IndexMut;

/// Two dimensional table
pub struct Table<T> {
    /// Width of the table.
    pub width: usize,
    /// Height of the table.
    pub height: usize,
    data: Vec<T>,
}

impl<T> Index<[usize; 2]> for Table<T> {
    type Output = T;

    /// Returns a reference to the value as (`x`, `y`).
    ///
    /// # Panics
    /// Panics if `x >= self.width` or `y >= self.height`.
    ///
    /// Panics if the `self` is not initialized
    fn index(&self, [x, y]: [usize; 2]) -> &T {
        if x < self.width && y < self.height {
            &self.data[self.width * y + x]
        } else {
            panic!(format!("{:?} is not a valid index", [x, y]))
        }
    }
}

impl<T> IndexMut<[usize; 2]> for Table<T> {
    /// Returns a mutable reference to the value as (`x`, `y`).
    ///
    /// # Panics
    /// Panics if `x >= self.width` or `y >= self.height`.
    ///
    /// Panics if the `self` is not initialized
    fn index_mut(&mut self, [x, y]: [usize; 2]) -> &mut T {
        if x < self.width && y < self.height {
            &mut self.data[self.width * y + x]
        } else {
            panic!(format!("{:?} is not a valid index", [x, y]))
        }
    }
}

impl<T> Table<T> {
    /// Allocates an empty table with `width` and `height`.
    ///
    /// The data inside is not initialized.
    pub fn new(width: usize, height: usize) -> Self {
        Table {
            width,
            height,
            data: Vec::with_capacity(width * height),
        }
    }
}

impl<T: Default> Table<T> {
    /// Initializes `self` with the default values of `T`.
    pub fn fill_default(&mut self) {
        self.data.resize_with(self.width * self.height, T::default);
    }
}

impl<T: Debug> Debug for Table<T> {
    /// Prints `self` as table, if `T`'s [Debug] implementation returns strings of fixed width.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for j in 0..self.height {
            for i in 0..self.width {
                write!(f, "{:?} ", self[[i, j]])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
