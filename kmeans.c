#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#define DIM 3

// Variáveis globais compartilhadas
pthread_mutex_t mutex;
pthread_cond_t cond;
int nthreads_global;

// Dados compartilhados
int *count_global;
double *sum_global;
int global_flips = 0;
int k_global;

int max_iteracoes = 1000; //garantir a convergência
int global_iteracao = 0;

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
            w->sum_local[c*DIM+0] = 0.0;
            w->sum_local[c*DIM+1] = 0.0;
            w->sum_local[c*DIM+2] = 0.0;
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

            w->count_local[cluster_global[i]]++;
            double *s = &w->sum_local[cluster_global[i] * DIM];
            s[0] += xi[0];
            s[1] += xi[1];
            s[2] += xi[2];
        }

        barreira(nthreads_global);

        if (tid == 0) {
            for (int c = 0; c < k_global; c++) {
                count_global[c] = 0;
                sum_global[c*DIM+0] = 0.0;
                sum_global[c*DIM+1] = 0.0;
                sum_global[c*DIM+2] = 0.0;
            }
            global_flips = 0;
            global_iteracao++;
        }

        barreira(nthreads_global);

        pthread_mutex_lock(&mutex);
        for (int c = 0; c < k_global; c++) {
            count_global[c] += w->count_local[c];
            double *sg = &sum_global[c * DIM];
            double *sl = &w->sum_local[c * DIM];
            sg[0] += sl[0];
            sg[1] += sl[1];
            sg[2] += sl[2];
        }
        global_flips += w->flips_local;
        pthread_mutex_unlock(&mutex);

        barreira(nthreads_global);

        if (tid == 0) {
            for (int c = 0; c < k_global; c++) {
                if (count_global[c] > 0) {
                    double inv = 1.0 / count_global[c];
                    mean_global[c*DIM + 0] = sum_global[c*DIM + 0] * inv;
                    mean_global[c*DIM + 1] = sum_global[c*DIM + 1] * inv;
                    mean_global[c*DIM + 2] = sum_global[c*DIM + 2] * inv;
                }
            }
        }

        barreira(nthreads_global);

        if (global_flips == 0 || global_iteracao >= max_iteracoes) break;
    }

    pthread_exit(NULL);
}

int main(int argc, char **argv) {
    int k, n;
    FILE *arquivo_entrada; // Ponteiro para o arquivo

    // --- 1. Checagem de argumentos (agora só 1: nome do arquivo) ---
    if (argc < 2) {
        fprintf(stderr, "Erro: informe o arquivo de entrada.\n");
        fprintf(stderr, "Uso: %s <arquivo_pontos>\n", argv[0]);
        return 1;
    }

    // --- 2. Pedir parâmetros interativamente ---
    printf("Digite o numero de threads: ");
    if (scanf("%d", &nthreads_global) != 1 || nthreads_global <= 0) {
        fprintf(stderr, "Numero de threads invalido\n");
        return 1;
    }

    printf("Digite o numero de clusters (k): ");
    if (scanf("%d", &k) != 1 || k <= 0) {
        fprintf(stderr, "Entrada invalida para k.\n");
        return 1;
    }
    k_global = k;

    // --- 3. Pedir N (número de pontos) interativamente ---
    printf("Digite o numero de pontos (n): ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Entrada invalida para n.\n");
        return 1;
    }
    n_global = n;

    // Validação
    if (k_global > n_global) {
        fprintf(stderr, "Erro: O numero de clusters (k=%d) nao pode ser maior que o numero de pontos (n=%d).\n", k_global, n_global);
        return 1;
    }

    // --- 4. Abrir o arquivo de pontos ---
    arquivo_entrada = fopen(argv[1], "r");
    if (arquivo_entrada == NULL) {
        perror("Erro ao abrir o arquivo de entrada");
        return 1;
    }

    // Alocações
    x_global = malloc(sizeof(double) * DIM * n);
    mean_global = malloc(sizeof(double) * DIM * k);
    cluster_global = malloc(sizeof(int) * n);

    if (!x_global || !mean_global || !cluster_global) {
        fprintf(stderr, "Erro de alocacao\n");
        fclose(arquivo_entrada);
        return 1;
    }

    // --- 5. Pedir centróides iniciais interativamente ---
    printf("\nDigite as coordenadas iniciais dos %d centroides (X Y Z):\n", k_global);
    for (int i = 0; i < k_global; i++) {
        printf("Centroide %d: ", i);
        if (scanf("%lf %lf %lf",
                  &mean_global[i*DIM],
                  &mean_global[i*DIM+1],
                  &mean_global[i*DIM+2]) != 3) {
            fprintf(stderr, "Erro na leitura das coordenadas do centroide.\n");
            fclose(arquivo_entrada); // Ainda precisa fechar o arquivo
            return 1;
        }
    }

    // --- 6. Ler pontos do arquivo ---
    // (Assume que o arquivo agora SÓ contém os pontos)
    printf("\nLendo %d pontos do arquivo '%s'...\n", n_global, argv[1]);
    for (int i = 0; i < n_global; i++) {
        if (fscanf(arquivo_entrada, "%lf %lf %lf",
                  &x_global[i*DIM],
                  &x_global[i*DIM+1],
                  &x_global[i*DIM+2]) != 3) {
            fprintf(stderr, "Erro ao ler dados do ponto %d no arquivo.\n", i);
            fclose(arquivo_entrada);
            return 1;
        }
        cluster_global[i] = 0; // Inicializa atribuição
    }

    // --- 7. Fechar o arquivo ---
    fclose(arquivo_entrada);
    
    // Alocações restantes
    count_global = malloc(sizeof(int) * k);
    sum_global = malloc(sizeof(double) * k * DIM);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);

    pthread_t *tid = malloc(sizeof(pthread_t) * nthreads_global);
    Worker *wk = malloc(sizeof(Worker) * nthreads_global);

    for (int t = 0; t < nthreads_global; t++) {
        wk[t].tid = t;
        wk[t].count_local = malloc(sizeof(int) * k);
        wk[t].sum_local = malloc(sizeof(double) * k * DIM);
    }

    clock_t inicio = clock();

    for (int t = 0; t < nthreads_global; t++) {
        pthread_create(&tid[t], NULL, worker_fn, &wk[t]);
    }

    for (int t = 0; t < nthreads_global; t++) {
        pthread_join(tid[t], NULL);
    }

    clock_t fim = clock();
    double tempo = (double)(fim - inicio) / CLOCKS_PER_SEC;

    // Saída dos centróides (stdout)
    printf("\nCentroides finais:\n");
    for (int i = 0; i < k; i++) {
        printf("Cluster %d => (%.2f, %.2f, %.2f)\n",
            i,
            mean_global[i*DIM], mean_global[i*DIM+1], mean_global[i*DIM+2]);
    }

    // Saída do tempo (stderr)
    fprintf(stderr, "\nTempo de execucao concorrente = %.6fs com %d threads\n",
           tempo, nthreads_global);

    // Limpeza de memória
    for (int t = 0; t < nthreads_global; t++) {
        free(wk[t].count_local);
        free(wk[t].sum_local);
    }
    free(wk);
    free(tid);
    free(count_global);
    free(sum_global);
    free(x_global);
    free(mean_global);
    free(cluster_global);
}