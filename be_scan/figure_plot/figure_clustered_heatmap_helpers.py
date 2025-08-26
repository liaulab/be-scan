
from be_scan.figure_plot.figure_classes import *

def kmeans_elbow(
    df,
    k_values,
    random_state=0,
    ):
    # ELBOW METHOD FOR K-MEANS ON RAW DATA #
    wcss = []  # Inertia
    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=random_state)
        kmeans.fit(df)  # Use raw data
        wcss.append(kmeans.inertia_)

    # PLOT ELBOW CURVE #
    plt.plot(k_values, wcss, marker='o')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
    plt.title('Elbow Method (Raw Data)')
    plt.xticks(k_values)
    plt.grid(True)
    plt.show()
    plt.close()

def kmeans_pca_scatterplot(
    df, k,
    random_state=0,
    ):

    # RUN K-MEANS ON RAW DATA #
    kmeans = KMeans(n_clusters=k, random_state=random_state)
    clusters = kmeans.fit_predict(df)

    # PCA FOR VISUALIZATION ONLY #
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(df)
    explained_variance = pca.explained_variance_ratio_
    print(f"Explained variance by PC1 and PC2: {explained_variance}")

    # Plot PCA projection colored by cluster
    fig, ax = plt.subplots()
    scatter = ax.scatter(X_pca[:, 0], X_pca[:, 1],
                         c=clusters, s=10, cmap='tab10')
    ax.set_xlabel(f"PC-1 ({explained_variance[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC-2 ({explained_variance[1]*100:.1f}% var)")
    plt.title("K-means Clustering (raw data), PCA for plotting")
    plt.colorbar(scatter, label='Cluster')
    plt.show()
    plt.close()

    return clusters

def hex_to_rgb(hex_color):
  hex_color = hex_color.lstrip('#')
  if len(hex_color) != 6:
      return None
  try:
      r = int(hex_color[0:2], 16)
      g = int(hex_color[2:4], 16)
      b = int(hex_color[4:6], 16)
      return (r/255, g/255, b/255)
  except ValueError:
      return None
