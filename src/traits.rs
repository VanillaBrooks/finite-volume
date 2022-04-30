use crate::prelude::*;

pub(crate) trait MethodMarker {
    type Method;

    fn create_method(n: usize, state: &State) -> Self::Method;
}

pub(crate) trait SolverMethod {
    fn solver_step(&mut self, q: &mut Matrix, f: &Matrix, state: &mut State);
}

pub(crate) trait SanityCheck{
    fn has_nans(&self)-> bool ;
    fn has_negs(&self)-> bool ;
    fn has_zero(&self)-> bool ;
}

impl <D>SanityCheck for ndarray::ArrayBase<ndarray::OwnedRepr<f64>, D> 
where D: ndarray::Dimension
{
    fn has_nans(&self) -> bool {
        self.iter().any(|x| x.is_nan())
    }
    fn has_negs(&self) -> bool {
        self.iter().any(|x| *x < 0.0)
    }
    fn has_zero(&self) -> bool {
        self.iter().any(|x| *x == 0.0)
    }
}

impl <D>SanityCheck for ndarray::ArrayBase<ndarray::ViewRepr<&f64>, D> 
where D: ndarray::Dimension
{
    fn has_nans(&self) -> bool {
        self.iter().any(|x| x.is_nan())
    }
    fn has_negs(&self) -> bool {
        self.iter().any(|x| *x < 0.0)
    }
    fn has_zero(&self) -> bool {
        self.iter().any(|x| *x == 0.0)
    }
}
