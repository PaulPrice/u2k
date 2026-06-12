from concurrent.futures import ProcessPoolExecutor


class MapCallable:
    """A wrapper around a callable that allows us to catch exceptions and
    re-raise them with additional context.
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, index, *args, **kwargs):
        try:
            return self.func(*args, **kwargs)
        except Exception as exc:
            raise RuntimeError(f"Error in map job with index {index}") from exc


class ProcessPool:
    """A convenience wrapper around `concurrent.futures.ProcessPoolExecutor`
    
    that allows us to wait for all jobs to finish when exiting a context
    manager.

    Parameters
    ----------
    workers : `int`
        The number of worker processes to use.
    """
    def __init__(self, workers: int):
        self.executor = ProcessPoolExecutor(workers)
        self.jobs = []
        self.results = []

    def clear(self):
        """Clear the list of jobs and results."""
        self.jobs = []
        self.results = []

    def submit(self, func, *args, **kwargs):
        """Submit a job to the process pool.

        You should use this method within a context manager.

        Parameters
        ----------
        func : `callable`
            The function to execute in the worker process.
        *args
            Positional arguments to pass to the function.
        **kwargs
            Keyword arguments to pass to the function.

        Returns
        -------
        job : `concurrent.futures.Future`
            A future representing the execution of the job.
        """
        job = self.executor.submit(func, *args, **kwargs)
        self.jobs.append(job)
        return job

    def map(self, func, *iterables, **kwargs):
        """Map a function over iterables in the process pool.

        You should NOT call this method within a context manager.

        Parameters
        ----------
        func : `callable`
            The function to execute in the worker processes.
        *iterables
            Iterables to map over. Each iterable should have the same length.
        **kwargs
            Keyword arguments to pass to the function.

        Returns
        -------
        results : list
            A list of results from applying the function to the iterables.
        """
        callable = MapCallable(func)
        with self:
            for ii, things in enumerate(zip(*iterables)):
                self.submit(callable, ii, *things, **kwargs)
        return self.results

    def __enter__(self):
        self.clear()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            for job in self.jobs:
                self.results.append(job.result())
        except Exception:
            for jj in self.jobs:
                jj.cancel()
            self.executor.shutdown(wait=True)
            self.executor = None
            raise
