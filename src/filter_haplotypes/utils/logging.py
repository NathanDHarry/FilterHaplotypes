"""
Logging utilities for FilterHaplotypes.
Sets up multi-level logging to console and file.
"""

import logging
import sys
import multiprocessing
from pathlib import Path
from logging.handlers import QueueHandler, QueueListener

def setup_logging(output_dir: Path):
    """
    Setup logging to both stdout (INFO) and log.txt (DEBUG) in output directory.
    Supports multiprocessing via a QueueListener.
    
    :param output_dir: Directory to save log.txt.
    :return: A Queue that can be passed to workers for logging.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "log.txt"
    
    # Formatter - strictly follow spec: timestamps, log levels, and module names
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Console handler (INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    # File handler (DEBUG)
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    
    # Queue for multiprocessing
    queue = multiprocessing.Manager().Queue(-1)
    
    # Listener in the main process
    listener = QueueListener(queue, console_handler, file_handler, respect_handler_level=True)
    listener.start()
    
    # Configure root logger in main process to use the queue
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    # Remove existing handlers
    for h in root.handlers[:]:
        root.removeHandler(h)
    root.addHandler(QueueHandler(queue))
    
    root.info(f"Logging initialized. Log file: {log_file}")
    
    return queue, listener

def worker_configurer(queue):
    """
    Configure a worker process to log to the central queue.
    """
    h = QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.DEBUG)
