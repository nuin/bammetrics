"""fifo_buffer(): the registry function checks its options, then raises UsageError.

FifoBuffer copies the standard input of a command to its standard output. A Python
call has no pipe to buffer and nothing to return, so it refuses and names the command.
"""
import pytest

import bammetrics

REFUSAL = "FifoBuffer copies the standard input of the command to its standard output"


def test_the_call_refuses_and_names_the_command():
    with pytest.raises(bammetrics.UsageError, match=REFUSAL) as raised:
        bammetrics.fifo_buffer()
    assert "bammetrics fifo-buffer" in str(raised.value)
    assert isinstance(raised.value, ValueError)
    with pytest.raises(bammetrics.UsageError, match=REFUSAL):
        bammetrics.fifo_buffer(buffer_size=1024, io_size=16, debug_frequency=1, name="lane 1",
                               verbosity="WARNING")


def test_options_are_checked_before_the_refusal():
    for keywords, message in [
        ({"buffer_size": 0}, "BUFFER_SIZE=0: a size from 1 to 2147483645 bytes is needed"),
        ({"buffer_size": "abc"}, "Problem constructing Integer from the string 'abc'"),
        ({"io_size": -1}, "IO_SIZE=-1: a size from 1 to 2147483645 bytes is needed"),
        ({"debug_frequency": 1.5}, "Problem constructing Integer"),
        ({"verbosity": "LOUD"}, "one of ERROR, WARNING, INFO, DEBUG"),
        ({"output": "copy.bin"}, "Unrecognized option: OUTPUT"),
        ({"input": "in.bin"}, "Unrecognized option: INPUT"),
        ({"stop_after": 1}, "STOP_AFTER"),
    ]:
        with pytest.raises(bammetrics.UsageError, match=message):
            bammetrics.fifo_buffer(**keywords)


def test_the_tool_is_listed_with_its_verb():
    tool = next(t for t in bammetrics.tools() if t == "FifoBuffer")
    assert tool.verb == "fifo-buffer"
    assert "UsageError" in tool.summary
