from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SyncManager
from multiprocessing import Manager
from contextlib import contextmanager
from dataclasses import dataclass


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


@dataclass
class Block:
    num_molecules: int
    ncores: int


def setup_blocks(ncores: int, num_molecules: int, mincores: int) -> list[Block]:
    blocks: list[Block] = []

    # Maximum and minimum number of parallel processes possible
    maxcores = ncores
    maxprocs = max(1, ncores // mincores)
    minprocs = 1

    # Distribute number of molecules among blocks
    # First (if possible) create the maximum number of parallel blocks (maxprocs) and distribute as many molecules as possible
    molecules_left = num_molecules
    if molecules_left >= maxprocs:
        p = maxprocs
        molecules_per_block = molecules_left // p
        for _ in range(p):
            blocks.append(Block(molecules_per_block, ncores // p))
        molecules_left -= molecules_per_block * p

    # While there are more than minprocs (1) molecules left find the optimal number of parallel blocks
    # Again distribute as many molecules per block as possible
    while molecules_left >= minprocs:
        p = max(
            [
                j
                for j in range(minprocs, maxprocs)
                if ncores % j == 0 and j <= molecules_left
            ]
        )
        molecules_per_block = molecules_left // p
        for _ in range(p):
            blocks.append(Block(molecules_per_block, ncores // p))
        molecules_left -= molecules_per_block * p

    # NOTE: using minprocs = 1 this is probably never true
    if molecules_left > 0:
        blocks.append(Block(molecules_left, maxcores))
        molecules_left -= molecules_left

    return blocks
