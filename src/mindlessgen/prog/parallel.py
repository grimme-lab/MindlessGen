from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SyncManager
from multiprocessing import Manager
from contextlib import contextmanager


@contextmanager
def setup_managers(max_workers: int, ncores: int):
    executor: ProcessPoolExecutor = ProcessPoolExecutor(max_workers=max_workers)
    manager: SyncManager = Manager()
    resource_manager: ResourceMonitor = ResourceMonitor(manager, ncores)
    try:
        yield executor, manager, resource_manager
    finally:
        executor.shutdown(False, cancel_futures=True)
        manager.shutdown()


class ResourceMonitor:
    def __init__(self, manager: SyncManager, ncores: int):
        self.__free_cores = manager.Value(int, ncores)
        self.__enough_cores = manager.Condition()

    @contextmanager
    def occupy_cores(self, ncores: int):
        try:
            with self.__enough_cores:
                self.__enough_cores.wait_for(lambda: self.__free_cores.value >= ncores)
                self.__free_cores.value -= ncores
            yield
        finally:
            with self.__enough_cores:
                self.__free_cores.value += ncores
                self.__enough_cores.notify()  # TODO: try this with notify_all instead
