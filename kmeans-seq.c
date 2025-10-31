#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIM 3
int max_iteracoes = 1000; //garantir a convergência
int global_iteracao = 0;

clock_t start, end;
double tempo_execucao;

// A função main agora aceita argumentos de linha de comando
int main(int argc, char *argv[]) {
    int i, j, k, n, c;
    double dmin, dx;
    double *x, *mean, *sum;
    int *cluster, *count, color;
    int flips;
    FILE *arquivo_pontos; // Ponteiro para o arquivo de entrada

    // --- 1. Checagem de Argumentos ---
    if (argc < 2) {
        fprintf(stderr, "Erro: informe o arquivo de entrada.\n");
        fprintf(stderr, "Uso: %s <arquivo_pontos>\n", argv[0]);
        return 1;
    }

    // --- 2. Leitura Interativa de K e N ---
    printf("Digite o numero de clusters (k): ");
    if (scanf("%d", &k) != 1 || k <= 0) {
        fprintf(stderr, "Entrada invalida para k.\n");
        return 1;
    }

    printf("Digite o numero de pontos (n): ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Entrada invalida para n.\n");
        return 1;
    }

    if (k > n) {
        fprintf(stderr, "Erro: k (%d) nao pode ser maior que n (%d).\n", k, n);
        return 1;
    }

    // --- 3. Abertura do Arquivo de Pontos ---
    arquivo_pontos = fopen(argv[1], "r");
    if (arquivo_pontos == NULL) {
        perror("Erro ao abrir o arquivo");
        return 1;
    }

    // --- 4. Alocação de Memória ---
    x = (double *)malloc(sizeof(double)*DIM*n);
    mean = (double *)malloc(sizeof(double)*DIM*k);
    sum= (double *)malloc(sizeof(double)*DIM*k);
    cluster = (int *)malloc(sizeof(int)*n);
    count = (int *)malloc(sizeof(int)*k);

    if (x == NULL || mean == NULL || sum == NULL || cluster == NULL || count == NULL) {
        fprintf(stderr, "Erro de alocacao de memoria.\n");
        fclose(arquivo_pontos);
        return 1;
    }

    // Inicializa clusters (como no original)
    for (i = 0; i<n; i++) 
        cluster[i] = 0;

    // --- 5. Leitura Interativa dos Centróides ---
    printf("\nDigite as coordenadas iniciais dos %d centroides (X Y Z):\n", k);
    for (i = 0; i<k; i++) {
        printf("Centroide %d: ", i);
        if (scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2) != 3) {
            fprintf(stderr, "Erro na leitura das coordenadas.\n");
            fclose(arquivo_pontos);
            return 1;
        }
    }

    // --- 6. Leitura dos Pontos (do Arquivo) ---
    printf("\nLendo %d pontos do arquivo '%s'...\n", n, argv[1]);
    for (i = 0; i<n; i++) {
        // Usa fscanf para ler do arquivo
        if (fscanf(arquivo_pontos, "%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2) != 3) {
            fprintf(stderr, "Erro ao ler o ponto %d do arquivo.\n", i);
            fprintf(stderr, "Certifique-se que o arquivo contem %d linhas de pontos.\n", n);
            fclose(arquivo_pontos);
            return 1;
        }
    }
    // Fecha o arquivo após a leitura
    fclose(arquivo_pontos);

    start = clock();
    // --- 7. Lógica K-Means (Inalterada) ---
    flips = n;
    while (flips>0 && global_iteracao < max_iteracoes) {
        flips = 0;
        for (j = 0; j < k; j++) {
            count[j] = 0; 
            for (i = 0; i < DIM; i++) 
                sum[j*DIM+i] = 0.0;
        }
        for (i = 0; i < n; i++) {
            dmin = -1; color = cluster[i];
            for (c = 0; c < k; c++) {
                dx = 0.0;
                for (j = 0; j < DIM; j++) 
                    dx +=  (x[i*DIM+j] - mean[c*DIM+j])*(x[i*DIM+j] - mean[c*DIM+j]);
                if (dx < dmin || dmin == -1) {
                    color = c;
                    dmin = dx;
                }
            }
            if (cluster[i] != color) {
                flips++;
                cluster[i] = color;
            }
        }

        for (i = 0; i < n; i++) {
            count[cluster[i]]++;
            for (j = 0; j < DIM; j++) 
                sum[cluster[i]*DIM+j] += x[i*DIM+j];
        }
        for (i = 0; i < k; i++) {
            // Adiciona checagem para evitar divisão por zero
            if (count[i] > 0) {
                for (j = 0; j < DIM; j++) {
                    mean[i*DIM+j] = sum[i*DIM+j]/count[i];
                }
            }
        }
        global_iteracao++;
    }
    end = clock();
    tempo_execucao = ((double)(end - start)) / CLOCKS_PER_SEC;

    // --- 8. Impressão dos Resultados (Inalterada) ---
    printf("\nCentroides finais:\n");
    for (i = 0; i < k; i++) {
        for (j = 0; j < DIM; j++)
            printf("%5.2f ", mean[i*DIM+j]);
        printf("\n");
    }

    printf("Tempo de execucao sequencial: %.6f segundos\n", tempo_execucao);

    #ifdef DEBUG
    for (i = 0; i < n; i++) {
        for (j = 0; j < DIM; j++)
            printf("%5.2f ", x[i*DIM+j]);
        printf("%d\n", cluster[i]);
    }
    #endif

    // --- 9. Limpeza da Memória ---
    free(x);
    free(mean);
    free(sum);
    free(cluster);
    free(count);
    
    return(0);
}