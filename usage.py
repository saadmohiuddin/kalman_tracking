from src import TrackGenerator

# Example showing how to generate tracks for Kalman Filter algorithm. Input for Saad
print("-" * 25, "Start", "-" * 25)

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
print("-" * 25, "Original Track", "-" * 25)
print(track)
print("-" * 50)

# reachability of original track - i.e how many layer does this track pass through
reachability = track.check_layer_reachability()
print(f"Reaches all Layer : {reachability['all_reachable']}")
print("-" * 50)
# print(reachability) # to reach all layers - all_reachable should be true.
# Saad - uncomment this if you want to see all info. But this is not necessary for your
# part of work
if not reachability["all_reachable"]:
    print("-" * 25, "Generating Track that reaches all layer", "-" * 25)
    # track with 100 points of X/Y/Z coordinates#
    hits = track.generate_track_with_all_hits(max_attempts=100)
    # for Aissata
    X, Y, Z, _, _ = track.evolve_track(
        time_steps=100
    )  # you can increase time_steps for granularity
    # just to show that it works
    print(X[:5], Y[:5], Z[:5], type(X), type(Y), type(Z))
    reachability = track.check_layer_reachability()

# this track hits all layers
radius = (33.25, 50.5, 88.5, 122.5)  # radius of pixel layer
hits = track.find_layer_intersections()
hits_layer1 = hits[radius[0]]
hits_layer2 = hits[radius[1]]
hits_layer3 = hits[radius[2]]
hits_layer4 = hits[radius[3]]
print("-" * 50)
print(f"Reaches all Layer : {reachability['all_reachable']}")
print("-" * 50)
print("XYZ coordinate of hit in layer 1", hits_layer1)  # in mm
print("XYZ coordinate of hit in layer 2", hits_layer2)
print("XYZ coordinate of hit in layer 3", hits_layer3)
print("XYZ coordinate of hit in layer 4", hits_layer4)

print("-" * 25, "End", "-" * 25)
