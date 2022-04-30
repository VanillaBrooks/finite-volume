pub use ndarray::s;
pub(crate) use crate::vector_calculations::recalculate_q;
pub(crate) use crate::vector_calculations::recalculate_f;
pub(crate) use crate::vector_calculations::recalculate_variables;
pub(crate) use crate::traits::*;
pub(crate) use crate::state::*;
pub(crate) use crate::solver_methods::*;

pub type Matrix = ndarray::Array2<f64>;
pub type Vector = ndarray::Array1<f64>;
