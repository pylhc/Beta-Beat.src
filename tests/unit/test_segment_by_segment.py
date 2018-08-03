import sys
import os
import shutil
import pytest
import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from segment_by_segment import segment_by_segment
from segment_by_segment.segment_by_segment import (
    SbsDefinitionError,
    Segment,
    SegmentModels,
    SegmentBeatings,
    GetLlmMeasurement,
)


CURRENT_DIR = os.path.dirname(__file__)


# _parse_elements() ###########################################################

def test_empty_parse_elements():
    with pytest.raises(SbsDefinitionError):
        segment_by_segment._parse_elements("")
    with pytest.raises(SbsDefinitionError):
        segment_by_segment._parse_elements("    ")


def test_correct_parse_elements():
    elem_str = "elem1,elem2"
    elements = segment_by_segment._parse_elements(elem_str)
    assert len(elements) == 2
    assert elements[0] == "elem1"
    assert elements[1] == "elem2"


def test_parse_elements_with_spaces():
    elem_str = "elem1, "
    elements = segment_by_segment._parse_elements(elem_str)
    assert len(elements) == 1
    assert "elem1" == elements[0]

    elem_str = "   elem1,  elem2  "
    elements = segment_by_segment._parse_elements(elem_str)
    assert len(elements) == 2
    assert "elem1" == elements[0]
    assert "elem2" == elements[1]


def test_parse_elements_with_duplicates():
    elem_str = "elem1,elem1"
    with pytest.raises(SbsDefinitionError):
        segment_by_segment._parse_elements(elem_str)


# _parse_segments() ###########################################################

def test_parse_segment_empty():
    segment_str = ""
    with pytest.raises(SbsDefinitionError):
        segment_by_segment._parse_segments(segment_str)

    segment_str = "     "
    with pytest.raises(SbsDefinitionError):
        segment_by_segment._parse_segments(segment_str)


def test_parse_correct_segment():
    segment_str = "name,start,end"
    segments = segment_by_segment._parse_segments(segment_str)
    assert len(segments) == 1
    assert segments[0].name == "name"
    assert segments[0].start == "start"
    assert segments[0].end == "end"

    segment_str = "name,start,end;"
    segments = segment_by_segment._parse_segments(segment_str)
    assert len(segments) == 1
    assert segments[0].name == "name"
    assert segments[0].start == "start"
    assert segments[0].end == "end"

    segment_str = "name,start,end;name2,start2,end2"
    segments = segment_by_segment._parse_segments(segment_str)
    assert len(segments) == 2
    assert segments[0].name == "name"
    assert segments[0].start == "start"
    assert segments[0].end == "end"
    assert segments[1].name == "name2"
    assert segments[1].start == "start2"
    assert segments[1].end == "end2"


# _there_are_duplicated_names() ###############################################

def test_there_are_duplicated_names():
    seg1 = Segment("name1", "", "")
    seg2 = Segment("name2", "", "")
    elements = ["elem1", "elem2"]
    assert not segment_by_segment._there_are_duplicated_names(
        [seg1, seg2], elements
    )
    elements = ["elem1", "name2"]
    assert segment_by_segment._there_are_duplicated_names(
        [seg1, seg2], elements
    )


def test_parse_segments_with_spaces():
    segment_str = "name,   start, end;   name2,start2,end2;  "
    segments = segment_by_segment._parse_segments(segment_str)
    assert len(segments) == 2
    assert segments[0].name == "name"
    assert segments[0].start == "start"
    assert segments[0].end == "end"
    assert segments[1].name == "name2"
    assert segments[1].start == "start2"
    assert segments[1].end == "end2"


# improve_segment() ###########################################################


def test_improve_segment_removes_forward(_test_df_just_names):
    seg = Segment("IP1", "BPM3", "BPM4")
    new_seg = segment_by_segment.improve_segment(
        seg,
        _test_df_just_names,
        _test_df_just_names,
        lambda name, _: name != "BPM4"
    )
    assert new_seg.name == seg.name
    assert new_seg.start == seg.start
    # BPM4 is wrong so it should jump it:
    assert new_seg.end == "BPM5"


def test_improve_segment_removes_backwards(_test_df_just_names):
    seg = Segment("IP1", "BPM3", "BPM4")
    new_seg = segment_by_segment.improve_segment(
        seg,
        _test_df_just_names,
        _test_df_just_names,
        lambda name, _: name != "BPM3"
    )
    assert new_seg.name == seg.name
    # BPM3 is wrong so it should jump it:
    assert new_seg.start == "BPM2"
    assert new_seg.end == seg.end


def test_improve_segment_removes_and_wraps(_test_df_just_names):
    seg = Segment("IP1", "BPM3", "BPM4")
    new_seg = segment_by_segment.improve_segment(
        seg,
        _test_df_just_names,
        _test_df_just_names,
        lambda name, _: name not in ("BPM1", "BPM2", "BPM3")
    )
    assert new_seg.name == seg.name
    # BPMs 1, 2 or 3 are wrong so it should wrap around to BPM6:
    assert new_seg.start == "BPM6"
    assert new_seg.end == seg.end


def test_improve_segment_raises_on_all_bad(_test_df_just_names):
    seg = Segment("IP1", "BPM3", "BPM4")
    with pytest.raises(SbsDefinitionError):
        segment_by_segment.improve_segment(
            seg,
            _test_df_just_names,
            _test_df_just_names,
            # All BPMs are bad:
            lambda name, _: False
        )


def test_improve_segment_same_if_all_good(_test_df_just_names):
    seg = Segment("IP1", "BPM3", "BPM4")
    new_seg = segment_by_segment.improve_segment(
        seg,
        _test_df_just_names,
        _test_df_just_names,
        # All BPMs are good:
        lambda name, _: True
    )
    assert new_seg.name == seg.name
    assert new_seg.start == seg.start
    assert new_seg.end == seg.end


def test_improve_segment_for_elements():
    seg = Segment("element", "element", "element")
    test_mdl = pd.DataFrame({"NAME": ["BPM1", "BPM2", "BPM3",
                                      "element",
                                      "BPM4", "BPM5", "BPM6"]})
    new_seg = segment_by_segment.improve_segment(
        seg,
        test_mdl,
        None,
        lambda name, _: name.startswith("BPM")
    )
    assert new_seg.name == seg.name
    assert new_seg.element == seg.element
    assert new_seg.start == "BPM3"
    assert new_seg.end == "BPM4"


# GetLlmMeasurement ###########################################################

def test_measurement_empty_dir(_meas_dir_empty):
    meas = GetLlmMeasurement(_meas_dir_empty)
    with pytest.raises(IOError):
        meas.beta_x
    with pytest.raises(IOError):
        meas.beta["x"]
    with pytest.raises(IOError):
        meas.coupling


def test_measurement_load_free(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    meas.beta["x"].BETX
    meas.beta_y.BETY


def test_measurement_load_free2(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    meas.phase_x.PHASEX
    meas.phase_y.PHASEY


def test_measurement_load_normal(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    meas.disp_x.DX
    meas.disp_y.DY


def test_measurement_load_no_plane(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    meas.coupling.F1001W


def test_measurement_load_norm_dips(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    meas.norm_disp.NDX


def test_measurement_writes_double_plane(_test_df_just_names, _test_dir):
    meas = GetLlmMeasurement(_test_dir)
    meas.allow_write = True
    meas.beta_x = (_test_df_just_names, "_free")
    meas.beta["y"] = (_test_df_just_names, "")
    assert os.path.isfile(
        os.path.join(CURRENT_DIR, "_test", "getbetax_free.out"))
    assert os.path.isfile(
        os.path.join(CURRENT_DIR, "_test", "getbetay.out"))


def test_measurement_doesnt_write_if_not_allowed(_test_df_just_names, _test_dir):
    meas = GetLlmMeasurement(_test_dir)
    meas.allow_write = False
    meas.beta_x = (_test_df_just_names, "_free")
    assert not os.path.isfile(
        os.path.join(CURRENT_DIR, "_test", "getbetax_free.out"))


def test_measurement_maybe_doesnt_call_on_error(_meas_dir_empty):
    meas = GetLlmMeasurement(_meas_dir_empty)
    def assert_false(args):
        assert False
    meas.maybe_call.beta_x(assert_false, "whatever")


def test_measurement_maybe_calls_on_success(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    def raise_to_be_catched(tfs_file):
        tfs_file.BETX  # Check that the columns can be accessed
        raise _KnownError()
    with pytest.raises(_KnownError):
        meas.maybe_call.beta_x(raise_to_be_catched)


def test_measurement_maybe_accessible(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    def raise_to_be_catched(tfs_file):
        tfs_file.BETY  # Check that the columns can be accessed
        raise _KnownError()
    with pytest.raises(_KnownError):
        meas.maybe_call.beta_y(raise_to_be_catched)
    with pytest.raises(_KnownError):
        meas.maybe_call.beta["y"](raise_to_be_catched)


def test_measurement_maybe_more_args(_meas_dir):
    meas = GetLlmMeasurement(_meas_dir)
    def test_funct(tfs_file, number, kwnumber=None):
        tfs_file.BETY  # Check that the columns can be accessed
        return number, kwnumber
    rnumber, rkwnumber = meas.maybe_call.beta["y"](test_funct, 5, 7)
    assert rnumber == 5
    assert rkwnumber == 7


# SegmentModels ##############################################################

def test_segment_models_empty_dir():
    seg = Segment("IP1", "", "")
    seg_mdls = SegmentModels(os.path.join(CURRENT_DIR, "meas_dir_empty"), seg)
    with pytest.raises(IOError):
        seg_mdls.front

def test_segment_models_returns_properly(_segment_models_dir):
    seg = Segment("IP1", "", "")
    seg_mdls = SegmentModels(_segment_models_dir, seg)
    seg_mdls.front.BETX
    seg_mdls.back.BETX
    seg_mdls.front_corrected.BETX
    seg_mdls.back_corrected.BETX


# SegmentBeatings #############################################################

def test_segment_beatings_writes_double_plane(_test_df_just_names, _test_dir):
    seg_beatings = SegmentBeatings(_test_dir, "test_seg")
    seg_beatings.allow_write = True
    seg_beatings.beta_phase_x = _test_df_just_names
    seg_beatings.beta_phase_y = _test_df_just_names
    assert os.path.isfile(
        os.path.join(CURRENT_DIR, "_test", "sbsbetabeatingx_test_seg.out"))
    assert os.path.isfile(
        os.path.join(CURRENT_DIR, "_test", "sbsbetabeatingy_test_seg.out"))


# Utilities ###################################################################

@pytest.fixture()
def _test_dir():
    test_dir = os.path.join(CURRENT_DIR, "_test")
    os.mkdir(test_dir)
    try:
        yield test_dir
    finally:
        shutil.rmtree(test_dir)


@pytest.fixture()
def _test_df_just_names():
    return pd.DataFrame({"NAME": ["BPM1", "BPM2", "BPM3",
                                  "BPM4", "BPM5", "BPM6"]})


@pytest.fixture()
def _segment_models_dir():
    return os.path.join(CURRENT_DIR, "..", "inputs", "sbs_models")


@pytest.fixture()
def _meas_dir():
    return os.path.join(CURRENT_DIR, "..", "inputs",
                        "getllm_results", "hllhc_sim")


@pytest.fixture()
def _meas_dir_empty():
    return os.path.join(CURRENT_DIR, "..", "inputs",
                        "getllm_results", "empty")


class _KnownError(Exception):
    pass
