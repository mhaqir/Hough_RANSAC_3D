# Written By Mohammad Haghir (mhaghir@iu.edu) unless cited otherwise

import numpy as np


# =========================== 
# ------- 2d to 3d ---------- 
# =========================== 
# https://github.com/mileyan/pseudo_lidar with modification
def project_image_to_rect(uv_depth, u, v, d):
    ''' Input: nx3 first two channels are uv, 3rd channel
               is depth in rect camera coord.
        Output: nx3 points in rect camera coord.
    '''
    # c_u, c_v, f_u, f_v, b_x, b_y = 609, 172, 721, 721, 44/(-721), 0.2/(-721)
    c_u, c_v, f_u, f_v, b_x, b_y = u, v, d, d, 44/(-721), 0.2/(-721)

    n = uv_depth.shape[0]
    x = ((uv_depth[:, 0] - c_u) * uv_depth[:, 2]) / f_u + b_x
    y = ((uv_depth[:, 1] - c_v) * uv_depth[:, 2]) / f_v + b_y
    pts_3d_rect = np.zeros((n, 3))
    pts_3d_rect[:, 0] = x
    pts_3d_rect[:, 1] = y
    pts_3d_rect[:, 2] = uv_depth[:, 2]
    return pts_3d_rect


def project_depth_to_points(depth, u, v, d): #, max_high
    rows, cols = depth.shape
    c, r = np.meshgrid(np.arange(cols), np.arange(rows))
    points = np.stack([c, r, depth])
    points = points.reshape((3, -1))
    points = points.T
    cloud = project_image_to_rect(points, u, v, d)
    # valid = (cloud[:, 0] >= 0) & (cloud[:, 2] < max_high)
    return cloud#[valid]


def animate(frames, anim_name):
  fig = plt.figure('fig')
  imgs = [[plt.imshow(f)] for f in frames]  # , cmap='gray'
  anim = animation.ArtistAnimation(fig, imgs, interval = 30, blit = True, repeat_delay=0)
  writegif = animation.PillowWriter(fps = 30)
  anim.save('{}.gif'.format(anim_name), writegif)



def animate3D(x, y, z, DIR):
  fig = plt.figure()
  ax = Axes3D(fig)
  def init():
    ax.plot_trisurf(x, y, z, cmap = 'viridis', edgecolor = 'none')
    return fig,
  def animate(i):
    ax.view_init(elev=40, azim=i)
    return fig,

  anim = animation.FuncAnimation(fig, animate, init_func=init,
                                 frames=360, interval=20, blit=True)
  anim.save('{}.mp4'.format(DIR[DIR.rfind('/') + 1:]),\
             fps=30, extra_args=['-vcodec', 'libx264'])

# p1 and p2 are two equal length array
# if the first three elements of at least
# one of the arrays are all zero, then this
# coparision doesn't work
def compare_planes(p1, p2):
  if np.sum(p1[0:3] != 0) == 0 or np.sum(p2[0:3] != 0) == 0:
    return None
  else: 
    nom = np.dot(p1[0:3], p2[0:3])
    denom = np.sqrt(np.sum(np.power(p1[0:3], 2)))*np.sqrt(np.sum(np.power(p2[0:3], 2)))#round(, 5)
    # print('nom/denom:', nom/denom)
    agl = np.arccos(nom/denom)
    if agl == 0:
      if p2[0] != 0:
        r = p1[0]/p2[0]
      elif p2[1] != 0:
        r = p1[1]/p2[1]
      else:
        r = p1[2]/p2[2]
      d = np.abs(p2[3] - p1[3]/r) / np.sqrt(np.sum(np.power(p2[0:3], 2)))
      return True, d
    #   print('d:', d)
    # print('nom:', nom)
    # print('denom:', denom)
    return False, agl * 180 / np.pi

