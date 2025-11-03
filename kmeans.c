#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#define DIM 3

pthread_mutex_t mutex;
pthread_cond_t cond;
int nthreads_global;

// Dados compartilhados
int *count_global;
double *sum_global;
int global_flips = 0;
int k_global;
int max_iteracoes = 1000;
int global_iteracao = 0;
volatile int convergiu = 0;

// Dados de entrada
double *x_global;
double *mean_global;
int *cluster_global;
int n_global;

typedef struct {
    int tid;
    int *count_local;
    double *sum_local;
    int flips_local;
} Worker;

Worker *wk_global;

void barreira(int nthreads) {
    static int bloqueadas = 0;
    pthread_mutex_lock(&mutex);
    if (bloqueadas == (nthreads - 1)) {
        bloqueadas = 0;
        pthread_cond_broadcast(&cond);
    } else {
        bloqueadas++;
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}

static inline double dist2_3d(const double *a, const double *b) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

void *worker_fn(void *arg) {
    Worker *w = (Worker *)arg;
    int tid = w->tid;

    int n_per_thread = (n_global + nthreads_global - 1) / nthreads_global;
    int start = tid * n_per_thread;
    int end = start + n_per_thread;
    if (end > n_global) end = n_global;

    for (;;) {
        for (int c = 0; c < k_global; c++) {
            w->count_local[c] = 0;
            for (int d = 0; d < DIM; d++)
                w->sum_local[c*DIM + d] = 0.0;
        }
        w->flips_local = 0;

        for (int i = start; i < end; i++) {
            const double *xi = &x_global[i * DIM];
            int old_c = cluster_global[i];
            double best = INFINITY;
            int best_c = 0;
            for (int c = 0; c < k_global; c++) {
                double d2 = dist2_3d(xi, &mean_global[c * DIM]);
                if (d2 < best) {
                    best = d2;
                    best_c = c;
                }
            }
            if (best_c != old_c) {
                cluster_global[i] = best_c;
                w->flips_local++;
            }
            w->count_local[best_c]++;
            for (int d = 0; d < DIM; d++)
                w->sum_local[best_c*DIM + d] += xi[d];
        }

        barreira(nthreads_global);

        int clusters_por_thread = (k_global + nthreads_global - 1) / nthreads_global;
        int c_start = tid * clusters_por_thread;
        int c_end = c_start + clusters_por_thread;
        if (c_end > k_global) c_end = k_global;

        for (int c = c_start; c < c_end; c++) {
            int count = 0;
            double sum[DIM] = {0};
            for (int t = 0; t < nthreads_global; t++) {
                count += wk_global[t].count_local[c];
                for (int d = 0; d < DIM; d++)
                    sum[d] += wk_global[t].sum_local[c*DIM + d];
            }
            count_global[c] = count;
            for (int d = 0; d < DIM; d++)
                sum_global[c*DIM + d] = sum[d];
        }

        barreira(nthreads_global);

        // Recalcular centróides
        for (int c = c_start; c < c_end; c++) {
            if (count_global[c] > 0) {
                double inv = 1.0 / count_global[c];
                for (int d = 0; d < DIM; d++)
                    mean_global[c*DIM + d] = sum_global[c*DIM + d] * inv;
            }
        }

        barreira(nthreads_global);

        if (tid == 0) {
            global_iteracao++;
            global_flips = 0;
            for (int t = 0; t < nthreads_global; t++)
                global_flips += wk_global[t].flips_local;
            convergiu = (global_flips == 0 || global_iteracao >= max_iteracoes);
        }

        barreira(nthreads_global);

        if (convergiu) break;
    }

    pthread_exit(NULL);
}

int main(int argc, char **argv) {
    int k, n;
    FILE *arquivo_entrada;

    if (argc < 2) {
        fprintf(stderr, "Uso: %s <arquivo_pontos>\n", argv[0]);
        return 1;
    }

    printf("Digite o numero de threads: ");
    scanf("%d", &nthreads_global);
    printf("Digite o numero de clusters (k): ");
    scanf("%d", &k);
    k_global = k;
    printf("Digite o numero de pontos (n): ");
    scanf("%d", &n);
    n_global = n;
    if (k_global > n_global) {
        fprintf(stderr, "Erro: k > n.\n");
        return 1;
    }

    arquivo_entrada = fopen(argv[1], "r");
    if (!arquivo_entrada) {
        perror("Erro ao abrir arquivo");
        return 1;
    }

    x_global = malloc(sizeof(double) * DIM * n);
    mean_global = malloc(sizeof(double) * DIM * k);
    cluster_global = malloc(sizeof(int) * n);

    printf("\nDigite as coordenadas iniciais dos %d centroides (X Y Z):\n", k_global);
    for (int i = 0; i < k_global; i++) {
        printf("Centroide %d: ", i);
        scanf("%lf %lf %lf",
              &mean_global[i*DIM],
              &mean_global[i*DIM+1],
              &mean_global[i*DIM+2]);
    }

    printf("\nLendo %d pontos de '%s'...\n", n_global, argv[1]);
    for (int i = 0; i < n_global; i++) {
        fscanf(arquivo_entrada, "%lf %lf %lf",
               &x_global[i*DIM],
               &x_global[i*DIM+1],
               &x_global[i*DIM+2]);
        cluster_global[i] = 0;
    }
    fclose(arquivo_entrada);

    count_global = malloc(sizeof(int) * k);
    sum_global = malloc(sizeof(double) * k * DIM);
    wk_global = malloc(sizeof(Worker) * nthreads_global);
    pthread_t *tid = malloc(sizeof(pthread_t) * nthreads_global);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);

    for (int t = 0; t < nthreads_global; t++) {
        wk_global[t].tid = t;
        wk_global[t].count_local = malloc(sizeof(int) * k);
        wk_global[t].sum_local = malloc(sizeof(double) * k * DIM);
    }

    clock_t inicio = clock();

    for (int t = 0; t < nthreads_global; t++)
        pthread_create(&tid[t], NULL, worker_fn, &wk_global[t]);

    for (int t = 0; t < nthreads_global; t++)
        pthread_join(tid[t], NULL);

    clock_t fim = clock();
    double tempo = (double)(fim - inicio) / CLOCKS_PER_SEC;

    if (global_iteracao >= max_iteracoes)
        printf("\nAviso: atingiu o limite de %d iteracoes.\n", max_iteracoes);

    printf("\nCentroides finais:\n");
    for (int i = 0; i < k_global; i++)
        printf("Cluster %d => (%.2f, %.2f, %.2f)\n",
               i,
               mean_global[i*DIM],
               mean_global[i*DIM+1],
               mean_global[i*DIM+2]);

    fprintf(stderr, "\nTempo total = %.6fs com %d threads\n",
           tempo, nthreads_global);

    for (int t = 0; t < nthreads_global; t++) {
        free(wk_global[t].count_local);
        free(wk_global[t].sum_local);
    }
    free(wk_global);
    free(tid);
    free(count_global);
    free(sum_global);
    free(x_global);
    free(mean_global);
    free(cluster_global);

    return 0;
}
