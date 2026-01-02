import numpy as np
from numpy.linalg import norm, cond
from numpy import array
#from simplex import SimplexNoise
from sine_wave_noise_3d import SineWaveNoise3D
#from numba import njit
from scipy.linalg import lstsq
from scipy.special import gamma
from scipy.spatial import cKDTree
from scipy.stats.qmc import Sobol
from math import log
from scipy.stats.qmc import PoissonDisk


def average_neighborhood_distance(n):
    """ Calculate the average distance between points in a neighborhood of size n."""
    m = 3**n-1  # Total number of points in the neighborhood
    d_avg = (1/m) * sum(2**i * binomial_coefficient(n, i) * i**0.5 for i in range(1,n+1))
    return d_avg

#@njit
def binomial_coefficient(n, m):
    if m < 0 or m > n:
        return 0
    if m == 0 or m == n:
        return 1
    m = min(m, n - m)  # Use symmetry
    c = 1
    for i in range(m):
        c = c * (n - i) // (i + 1)
    return c

def compute_rdf(points, max_radius, num_bins):
    """ Compute the radial distribution function (RDF) for a point cloud in n-dimensional space.
    Parameters:
    points (numpy.ndarray): The point cloud, an array of shape (N, D) where N is the number of points and D is the dimensionality.
    max_radius (float): The maximum radius to consider for the RDF.
    num_bins (int): The number of bins to use for the RDF.
    Returns:
    tuple: A tuple containing the RDF values and the corresponding radii."""
    N, D = points.shape
    radii = np.linspace(0, max_radius, num_bins + 1)
    rdf_values = np.zeros(num_bins)

    print("building KDTree", flush=True)
    tree = cKDTree(points)
    print("Finding pairs in KDTree", flush=True)
    pairs = tree.query_pairs(r=max_radius, output_type='ndarray')
    print("distances", flush=True)
    # Use vectorized computation for distances
    distances = np.linalg.norm(points[pairs[:,0]] - points[pairs[:,1]], axis=1)
    print("distances len: ", len(distances), flush=True)
    print("distances histogram", flush=True)
    hist, _ = np.histogram(distances, bins=radii)
    print("min distance", flush=True)
    min_dist = np.min(distances)

    # Correct shell volumes with D-ball prefactor
    prefactor = np.pi**(D/2) / gamma(D/2 + 1)
    shell_volumes = prefactor * np.diff(radii**D)
    density = N / (1.0 ** D) # Assuming unit hypercube
    rdf_values = hist / (shell_volumes * density * N)

    print("finding median and mean distancs", flush=True)
    # Query k=2 because the closest point to each point is itself (distance 0)
    dists, idxs = tree.query(points, k=2)
    # dists[:, 1] is the distance to the nearest neighbor (excluding self)
    avg_dist = np.mean(dists[:, 1])
    med_dist = np.median(dists[:, 1])
    print("done", flush=True)

    return rdf_values, radii[:-1], min_dist, avg_dist, med_dist


#@njit
def nD_cross_product(*vectors: list[np.ndarray]) -> np.ndarray:
    """Compute the n-dimensional cross product of n-1 vectors in n-dimensional space.
    Parameters:
    vectors: n-1 nD vectors passed as numpy arrays of type float64.
    Returns:
    numpy.ndarray: The resulting n-dimensional vector of shape (n,)."""
    n = len(vectors[0])
    if len(vectors) != n - 1:
        raise ValueError("Number of vectors must be one less than their dimension.")
    # Initialize the output vector
    output = np.zeros(n, dtype=np.float64)
    # Compute the sub-determinants for each coordinate
    for i in range(n):
        sub_matrix = np.zeros((n - 1, n - 1))
        for j,vec in enumerate(vectors):
            sub_matrix[j, :] = np.delete(vec, i)
        output[i] = (-1)**i * np.linalg.det(sub_matrix)
    return output


class DFVN_trace:
    """Class to perform n-dimensional divergence-free vector field tracing based on Simplex Noise."""
    def __init__(self, seed=42, dimensions=3, scale=1.0):
        """Initialize the DFVN trace with a given seed, dimensions, and scale."""
        self.seed = seed
        self.dimensions = dimensions
        self.noise = SineWaveNoise3D(seed=seed, num_waves=64, frequency_scale=1.0)
        self.scale = scale
        np.random.seed(seed)
        self.offsets = [np.random.uniform(0, 100.0 * self.scale, self.dimensions) for _ in range(self.dimensions - 1)]
 
    def curl_noise(self, p):
        scaling = [ 1.0, 1.1, 1.23, 2.16, 8.57, 83.99, 2847.73, 728677.50, 305579156.49]
        s = scaling[self.dimensions-3]
        vecs = []
        nvals = []
        for d in range(self.dimensions - 1):
            n, v = self.noise.noise_and_gradient(*(p*self.scale + self.offsets[d]))
            nvals.append(n)
            vecs.append(np.array(v)*self.scale)
        return s*nD_cross_product(*vecs), nvals

    def project(self, p, _iso_vals, step_size):
        '''This function projects point p onto the intersections of the iso-contours of n-1 noise functions.
        We are thus exploiting that the DFVN vector field integral curves are precisely these intersections.
        The method is to use the gradients to construct a linear system whose solutions are precisely the
        lines of intersection of the linear approximations to the iso-contours.'''
        iso_vals = np.array(_iso_vals, dtype=np.float64)
        def linear_system(p):
            # Build linear system:
            # rows of A are the gradients of the noise
            # b is the target noise values minus current value.
            A = np.zeros((self.dimensions-1,self.dimensions), dtype=np.float64)
            b = iso_vals.copy()
            for d in range(self.dimensions - 1):
                n, g = self.noise.noise_and_gradient(*(p*self.scale + self.offsets[d]))
                A[d, :] = np.array(g) * self.scale
                b[d] -= n
            return A,b
        
        A,b = linear_system(p)
        for i in range(1,11): # Max 10 iterations -- a high number to ensure convergence
            # Find the displacement as the least norm solution
            # cond is the threshold on singular values.
            if cond(A) > 1e3:
                break
            disp = lstsq(A, b, cond=0.1)[0]
            d = norm(disp)
            if d > 0.1*step_size or d < 1e-10*step_size: # If we are very close to the target, stop
                break
            p = p + disp # obtain new position
            A,b = linear_system(p)
        return p
    


    # def project(self, p, nvals):
    #     '''This is the old project function which takes one noise function
    #     at a time. '''
    #     p_old = p.copy()
    #     for d in range(self.dimensions - 1):
    #         n, v = self.noise.noise_grad(p*self.scale + self.offsets[d])
    #         v = np.array(v) * self.scale
    #         p = p - ((n-nvals[d])/ (v@v))*v
    #     return p, 1, 1


    def dfvn_trace(self, p, _t=1.0):
        """Perform the DFVN trace on a point p with a given time step _t."""
        t = _t / self.scale
        c_vec, nvals = self.curl_noise(p)
        a = t * c_vec
        c_vec, _ = self.curl_noise(p + a/2)
        b = t * c_vec
        c_vec, _ = self.curl_noise(p + b/2)
        c = t * c_vec
        c_vec, _ = self.curl_noise(p + c)
        d = t * c_vec
        p += (a + 2*b + 2*c + d) / 6
        # return p, 0, 0.0
        p, iters, b_norm = self.project(p, nvals, t)
        return p, iters, b_norm

    def dfvn_multi_trace(self, p, _t=1.0, N=1, w_project=True):
        """Perform the DFVN trace on a point p with a given time step _t."""
        t = _t / self.scale
        p = np.array(p)
        _, nvals = self.curl_noise(p)
        path = [p]
        for _ in range(N):
            c_vec, _ = self.curl_noise(p)
            # a = t * c_vec
            # c_vec, _ = self.curl_noise(p + a/2)
            # b = t * c_vec
            # c_vec, _ = self.curl_noise(p + b/2)
            # c = t * c_vec
            # c_vec, _ = self.curl_noise(p + c)
            # d = t * c_vec
            # p = p + (a + 2*b + 2*c + d) / 6
            p = p + t * c_vec
            if w_project:
                p = self.project(p, nvals,t)
            path.append(p)
        return path
    
    def jitter(self, p, delta):
        """Apply jitter to a point p based on the number of points,k, along each dimension."""
        return p + np.random.uniform(-delta/2, delta/2, self.dimensions)
    
def sobol(k, n, delta):
    num = 2**int(round(log(k**n)/log(2)))
    sobol_seq = Sobol(d=n)
    points = sobol_seq.random(num)
    u_bounds=np.ones(n)*delta*k
    return points*u_bounds

def pds(k, n, delta):
    l_bounds=np.zeros(n)
    u_bounds=np.ones(n)*delta*k
    pds_seq = PoissonDisk(d=n, radius=delta, l_bounds=l_bounds, u_bounds=u_bounds)
    points = pds_seq.random(k**n)
    return points

def create_nd_grid(k, n, delta=None):
    if delta is None:
        delta=1/k
    grid = np.array(np.meshgrid(*[np.linspace(delta/2, 1-delta/2, k) for _ in range(n)])).T.reshape(-1, n)
    return grid

def filter_grid(grid):
    # Filter grid to remove points outside the unit cube
    filtered_grid = grid[np.all((grid >= 0) & (grid <= 1), axis=1)]
    return filtered_grid
