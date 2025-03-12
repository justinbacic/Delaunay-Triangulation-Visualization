class triangulation:
    def __init__(self):
        self.points = set()
    def insertPoint(self, point):
        self.points.add(point)
    def print(self):
        for i in self.points:
            print(i)
tri = triangulation()
tri.insertPoint((3,4))
tri.print()
tri.insertPoint((4,4))
tri.print()
