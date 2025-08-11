import json
import matplotlib.pyplot as plt

min_xs = []
min_ys = []
max_xs = []
max_ys = []
xs = []
ys = []
zs = []
with open('250m_mesh_2024_14.geojson', 'r', encoding='shift-jis') as f:
    data = json.load(f)
    for feature in data['features']:
        min_x = 99999
        min_y = 99999
        max_x = -99999
        max_y = -99999
        for c in feature['geometry']['coordinates'][0]:
            min_x = min(min_x, c[0])
            min_y = min(min_y, c[1])
            max_x = max(max_x, c[0])
            max_y = max(max_y, c[1])

        z = feature['properties']['PT00_2025']

        min_xs.append(min_x)
        min_ys.append(min_y)
        max_xs.append(max_x)
        max_ys.append(max_y)
        xs.append(min_x)
        xs.append(max_x)
        ys.append(min_y)
        ys.append(max_y)
        zs.append(z)

mx = max(xs), min(xs)
my = max(ys), min(ys)
mz = max(zs), min(zs)
rx = mx[0] - mx[1]
ry = my[0] - my[1]
rz = mz[0] - mz[1]
xs0 = []
ys0 = []
xs1 = []
ys1 = []
zs0 = []
for i in range(len(min_xs)):
    mix = (min_xs[i] - mx[1]) / rx * 10000
    miy = (min_ys[i] - my[1]) / ry * 10000
    max = (max_xs[i] - mx[1]) / rx * 10000
    may = (max_ys[i] - my[1]) / ry * 10000
    z = (zs[i] - mz[1]) / rz
    if z < 0.1:
        continue
    z += 0.5
    z *= 3
    w = (max - mix) * z
    h = (may - miy) * z
    cx = (mix + max) / 2
    cy = (miy + may) / 2
    x0 = cx - w / 2
    y0 = cy - h / 2
    x1 = cx + w / 2
    y1 = cy + h / 2
    xs0.append(x0)
    ys0.append(y0)
    xs1.append(x1)
    ys1.append(y1)
    zs0.append(zs[i])

with open('sample2daabb_0.csv', 'w') as f:
    for mix, miy, max, may, z in zip(xs0, ys0, xs1, ys1, zs0):
        f.write(f"{mix},{miy},{max},{may},{z}\n")


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(xs, ys, zs, s=0.1)

plt.scatter(xs0, ys0, s=0.5)#, s=zs)
plt.show()

