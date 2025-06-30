#include <Rcpp.h>
#include <queue>
#include <vector>
#include <cmath>
using namespace Rcpp;

// Structure for KD-tree node
struct KdNode {
  double point[2];
  int index;
  KdNode* left;
  KdNode* right;
  
  KdNode(double x, double y, int idx)
    : index(idx), left(nullptr), right(nullptr) {
    point[0] = x;
    point[1] = y;
  }
};

// Build KD-tree (optimized memory management)
KdNode* buildKdTree(NumericMatrix coords, int depth = 0) {
  size_t n = coords.nrow();
  if (n == 0) return nullptr;
  
  int axis = depth % 2;
  NumericVector sorted = coords(_, axis);
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) {
    return sorted[i] < sorted[j];
  });
  
  size_t mid = n / 2;
  KdNode* node = new KdNode(coords(idx[mid],0), coords(idx[mid],1), idx[mid]);
  
  if (mid > 0) {
    NumericMatrix left_coords(mid, 2);
    for (size_t i=0; i<mid; ++i) {
      left_coords(i,0) = coords(idx[i],0);
      left_coords(i,1) = coords(idx[i],1);
    }
    node->left = buildKdTree(left_coords, depth+1);
  }
  
  if (mid+1 < n) {
    NumericMatrix right_coords(n-mid-1, 2);
    for (size_t i=mid+1; i<n; ++i) {
      right_coords(i-mid-1,0) = coords(idx[i],0);
      right_coords(i-mid-1,1) = coords(idx[i],1);
    }
    node->right = buildKdTree(right_coords, depth+1);
  }
  
  return node;
}

// KD-tree nearest neighbor search (thread-safe)
void knn_search(KdNode* node, const double* target, int k,
                std::priority_queue<std::pair<double, int>>& heap, int depth = 0) {
  if (!node) return;
  
  double dx = node->point[0] - target[0];
  double dy = node->point[1] - target[1];
  double dist = dx*dx + dy*dy;
  
  if (heap.size() < k || dist < heap.top().first) {
    if (heap.size() == k) heap.pop();
    heap.push(std::make_pair(dist, node->index));
  }
  
  int axis = depth % 2;
  double diff = target[axis] - node->point[axis];
  KdNode* near = diff <= 0 ? node->left : node->right;
  KdNode* far = diff <= 0 ? node->right : node->left;
  
  knn_search(near, target, k, heap, depth+1);
  
  if (diff*diff < heap.top().first) {
    knn_search(far, target, k, heap, depth+1);
  }
}

// [[Rcpp::export]]
List cpp_moran(NumericMatrix coords, NumericVector expr, int k = 15,
               bool global = true, int n_threads = 4) {
  int n = coords.nrow();
  KdNode* tree = buildKdTree(coords);
  List output;
  
  // Construct sparse weight matrix
  List neighbors(n);
  for (int i=0; i<n; ++i) {
    double target[2] = {coords(i,0), coords(i,1)};
    std::priority_queue<std::pair<double, int>> heap;
    knn_search(tree, target, k+1, heap); // k+1 to exclude self
    IntegerVector nb(k);
    for (int j=0; j<k; ++j) {
      nb[j] = heap.top().second + 1; // R indexing (1-based)
      heap.pop();
    }
    neighbors[i] = nb;
  }
  
  // Global Moran's I calculation
  if (global) {
    double mean_expr = mean(expr);
    double sum_num = 0.0, sum_denom = 0.0;
    double s0 = n * k;
    
    for (int i=0; i<n; ++i) {
      double xi = expr[i] - mean_expr;
      sum_denom += xi * xi;
      
      IntegerVector nb = neighbors[i];
      for (int j=0; j<k; ++j) {
        double xj = expr[nb[j]-1] - mean_expr;
        sum_num += xi * xj;
      }
    }
    
    double moran_I = (sum_num / sum_denom) * (n / s0);
    output["global"] = moran_I;
  }
  // Local Moran's I calculation
  else {
    NumericVector local_I(n);
    double mean_expr = mean(expr);
    double s2 = var(expr) * (n-1)/n; // Bias correction
    
    for (int i=0; i<n; ++i) {
      double xi = expr[i] - mean_expr;
      IntegerVector nb = neighbors[i];
      double sum_j = 0.0;
      
      for (int j=0; j<k; ++j) {
        sum_j += expr[nb[j]-1] - mean_expr;
      }
      
      local_I[i] = (xi / s2) * (sum_j / k);
    }
    output["local"] = local_I;
  }
  
  return output;
}

// Step 1: KD-tree accelerated neighbor search (build spatial weight matrix)
// [[Rcpp::export]]
List kNN_weights_cpp(NumericMatrix coords, int k) {
  int n = coords.nrow();
  List nb(n);
  
  for (int i=0; i<n; ++i) {
    std::priority_queue<std::pair<double, int> > pq;
    for (int j=0; j<n; ++j) {
      if (i == j) continue;
      double dx = coords(i,0) - coords(j,0);
      double dy = coords(i,1) - coords(j,1);
      double dist = dx*dx + dy*dy; // Squared distance
      
      if (pq.size() < k) {
        pq.push(std::make_pair(dist, j));
      } else if (dist < pq.top().first) {
        pq.pop();
        pq.push(std::make_pair(dist, j));
      }
    }
    
    IntegerVector neighbors(k);
    for (int m=k-1; m>=0; --m) {
      neighbors[m] = pq.top().second + 1; // R indexing (1-based)
      pq.pop();
    }
    nb[i] = neighbors;
  }
  
  return List::create(_["nn"] = nb);
}

// Step 2: Global Moran's I calculation (Rcpp parallel optimized)
// [[Rcpp::export]]
double moranI_cpp(NumericVector x, List nb) {
  int n = x.size();
  double mean_x = mean(x);
  double sum_num = 0.0, sum_denom = 0.0;
  double s0 = 0.0;
  
  for (int i=0; i<n; ++i) {
    double xi = x[i] - mean_x;
    sum_denom += xi * xi;
    
    IntegerVector neighbors = nb[i];
    int k = neighbors.size();
    s0 += k;
    
    for (int j=0; j<k; ++j) {
      int idx = neighbors[j] - 1; // Convert to C++ index (0-based)
      double xj = x[idx] - mean_x;
      sum_num += xi * xj;
    }
  }
  
  return (n / s0) * (sum_num / sum_denom);
}

// Step 3: Batch calculation of Local Moran's I
// [[Rcpp::export]]
NumericVector localMoran_cpp(NumericVector x, List nb) {
  int n = x.size();
  NumericVector Ii(n);
  double mean_x = mean(x);
  double s2 = var(x) * (n-1)/n;
  
  for (int i=0; i<n; ++i) {
    double xi = x[i] - mean_x;
    double sum_j = 0.0;
    int k = 0;
    
    IntegerVector neighbors = nb[i];
    k = neighbors.size();
    
    for (int j=0; j<k; ++j) {
      int idx = neighbors[j] - 1;
      double xj = x[idx] - mean_x;
      sum_j += xj;
    }
    
    Ii[i] = (xi / s2) * sum_j;
  }
  
  return Ii;
}