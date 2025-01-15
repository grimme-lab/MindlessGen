from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SyncManager
from multiprocessing import Manager
from typing import Any


class ParallelManager:
    def __init__(self, max_workers: int, ncores: int):
        self.executor: ProcessPoolExecutor = ProcessPoolExecutor(
            max_workers=max_workers
        )
        self.manager: SyncManager = Manager()

        self.__free_cores = self.manager.Value(int, ncores)
        self.__enough_cores = self.manager.Condition()

    def __enter__(self):
        return self

    def __exit__(self, exc_type: Exception, exc_val: Any, exc_tb: Any):
        self.executor.shutdown(False, cancel_futures=True)
        self.manager.shutdown()

        return False  # Don't supress any exceptions

    def shutdown(self):
        self.executor.shutdown(False, cancel_futures=True)
        self.manager.shutdown()

    def occupy_cores(self, ncores: int):
        with self.__enough_cores:
            self.__enough_cores.wait_for(lambda: self.__free_cores.value >= ncores)
            self.__free_cores.value -= ncores

    def free_cores(self, ncores: int):
        with self.__enough_cores:
            self.__free_cores.value += ncores
            self.__enough_cores.notify()  # TODO: try this with notify_all instead
