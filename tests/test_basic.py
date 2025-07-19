"""Basic tests"""

import pytest


def test_basic_imports():
    """Test that basic imports work."""
    from src import TrackGenerator

    assert TrackGenerator is not None


def test_track_generator_creation():
    """Test that TrackGenerator is properly initialized"""
    from src import TrackGenerator

    initial_guess = {
        "x": 2,
        "y": 1,
        "z": 0,
        "px": 30,
        "py": 40,
        "pz": 60,
        "charge": -1,
    }

    track = TrackGenerator.TrackGenerator(**initial_guess)
    assert track is not None
    assert track.x == 2
    assert track.y == 1
    assert track.z == 0
    assert track.px == 30
    assert track.py == 40
    assert track.pz == 60
    assert track.charge == -1


if __name__ == "__main__":
    pytest.main([__file__])
