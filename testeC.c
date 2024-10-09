#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#define MAX 3

typedef struct
{
    bool sisPossivel;
    bool sisDeterminado;
    float solucao[MAX];
} ResultadoSistema;

typedef struct
{
    double lambda1;
    double lambda2;
} Autovalores;

typedef struct
{
    float vetor1[2];
    float vetor2[2];
} Autovetores;

// Função para calcular o determinante de uma matriz 2x2
float determinante2x2(float a[2][2])
{
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

float determinante3x3(float a[MAX][MAX])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
           a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

// Função para calcular autovalores de uma matriz 2x2
Autovalores calcular_autovalores_2x2(float A[2][2])
{
    Autovalores autovalores;
    float trace = A[0][0] + A[1][1];                      // Traço da matriz
    float determinant = determinante2x2(A);               // Determinante
    float discriminant = trace * trace - 4 * determinant; // Delta

    if (discriminant < 0)
    {
        printf("O discriminante é negativo, autovalores complexos não são suportados.\n");
        autovalores.lambda1 = 0;
        autovalores.lambda2 = 0;
    }
    else
    {
        autovalores.lambda1 = (trace + sqrt(discriminant)) / 2;
        autovalores.lambda2 = (trace - sqrt(discriminant)) / 2;
    }

    return autovalores;
}

// Função para calcular autovetores de uma matriz 2x2
Autovetores calcular_autovetores_2x2(float A[2][2])
{
    Autovetores autovetores;
    Autovalores autovalores = calcular_autovalores_2x2(A);

    // Para lambda1
    float lambda1 = autovalores.lambda1;
    float mat1[2][2] = {
        {A[0][0] - lambda1, A[0][1]},
        {A[1][0], A[1][1] - lambda1}};

    if (mat1[0][1] != 0)
    {
        autovetores.vetor1[0] = 1;
        autovetores.vetor1[1] = -mat1[0][0] / mat1[0][1];
    }
    else if (mat1[1][0] != 0)
    {
        autovetores.vetor1[0] = -mat1[1][1] / mat1[1][0];
        autovetores.vetor1[1] = 1;
    }
    else
    {
        autovetores.vetor1[0] = 1;
        autovetores.vetor1[1] = 0;
    }

    // Para lambda2
    float lambda2 = autovalores.lambda2;
    float mat2[2][2] = {
        {A[0][0] - lambda2, A[0][1]},
        {A[1][0], A[1][1] - lambda2}};

    if (mat2[0][1] != 0)
    {
        autovetores.vetor2[0] = 1;
        autovetores.vetor2[1] = -mat2[0][0] / mat2[0][1];
    }
    else if (mat2[1][0] != 0)
    {
        autovetores.vetor2[0] = -mat2[1][1] / mat2[1][0];
        autovetores.vetor2[1] = 1;
    }

    else
    {
        autovetores.vetor2[0] = 1;
        autovetores.vetor2[1] = 0;
    }

    return autovetores;
}

// Função para adicionar histórico
void adicionarHistorico(const char *formato, ...)
{
    FILE *arquivo = fopen("historico.txt", "a");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir o arquivo. \n");
        return;
    }

    va_list args;
    va_start(args, formato);
    vfprintf(arquivo, formato, args);
    va_end(args);

    fclose(arquivo);
}

// Função para imprimir o histórico
void imprimirHistorico()
{
    FILE *arquivo = fopen("historico.txt", "r");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir o arquivo de historico.\n");
        return;
    }

    char linha[256];
    printf("Historico:\n");
    while (fgets(linha, sizeof(linha), arquivo))
    {
        printf("%s", linha);
    }

    fclose(arquivo);
}

// Função para limpar o histórico
void limparHistorico()
{
    FILE *arquivo = fopen("historico.txt", "w");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir o arquivo para limpeza.\n");
        return;
    }
    fclose(arquivo);
    printf("Historico limpo com sucesso.\n");
}

// Função para imprimir a matriz
void imprimirMatriz(float a[MAX][MAX + 1], int lin)
{
    for (int i = 0; i < lin; i++)
    {
        for (int j = 0; j <= lin; j++)
        {
            printf("%7.2f ", a[i][j]);
            adicionarHistorico("%7.2f ", a[i][j]);
        }
        adicionarHistorico("\n");
        printf("\n");
    }
}

// Função para aplicar a eliminação de Gauss
ResultadoSistema eliminacaoDeGauss(float a[MAX][MAX + 1], int lin)
{
    ResultadoSistema resultado;
    resultado.sisPossivel = true;
    resultado.sisDeterminado = true;

    for (int i = 0; i < lin - 1; i++)
    {
        if (a[i][i] == 0)
        {
            bool trocaEfetuada = false;

            for (int k = i + 1; k < lin; k++)
            {
                if (a[k][i] != 0)
                {
                    for (int j = 0; j <= lin; j++)
                    {
                        float temp = a[i][j];
                        a[i][j] = a[k][j];
                        a[k][j] = temp;
                    }
                    trocaEfetuada = true;
                    break;
                }
            }

            if (!trocaEfetuada)
            {
                resultado.sisPossivel = false;
                resultado.sisDeterminado = false;
                printf("Erro: Sistema singular ou indeterminado.\n");
                return resultado;
            }
        }

        // Escalonamento
        for (int k = i + 1; k < lin; k++)
        {
            float fator = a[k][i] / a[i][i];
            for (int j = 0; j <= lin; j++)
            {
                a[k][j] = a[k][j] - fator * a[i][j];
            }
        }
    }
    return resultado;
}

// Função para resolver o sistema após a matriz estar escalonada
ResultadoSistema resolucaoSistema(float a[MAX][MAX + 1], int lin)
{
    ResultadoSistema resultado;
    resultado.sisPossivel = true;
    resultado.sisDeterminado = true;

    float x[MAX];
    for (int i = lin - 1; i >= 0; i--)
    {
        if (a[i][i] == 0)
        {
            if (a[i][lin] != 0)
            {
                resultado.sisPossivel = false;
                printf("Erro: Sistema impossível.\n");
                return resultado;
            }
            else
            {
                resultado.sisDeterminado = false;
                printf("Erro: Sistema indeterminado.\n");
                return resultado;
            }
        }

        x[i] = a[i][lin];
        for (int j = i + 1; j < lin; j++)
        {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }

    for (int i = 0; i < lin; i++)
    {
        resultado.solucao[i] = x[i];
    }

    printf("Solucao do sistema:\n");
    for (int i = 0; i < lin; i++)
    {
        printf("x%d = %7.2f\n", i + 1, x[i]);
        adicionarHistorico("x%d = %7.2f\n", i + 1, x[i]);
    }

    return resultado;
}

// Função para ler a matriz do usuário
void scanearMatriz(float a[MAX][MAX + 1], int lin)
{
    printf("Digite os coeficientes das equacoes (ex: 1x + 2y - 3z = 4) = (1 2 -3 4)\n");
    for (int i = 0; i < lin; i++)
    {
        for (int j = 0; j <= lin; j++)
        {
            scanf("%f", &a[i][j]);
        }
    }
}

// Função para verificar se os vetores formam uma base
void verificarBase(float vetores[MAX][MAX], int dim)
{
    float det;

    if (dim == 2)
    {
        det = determinante2x2(vetores);
        if (det != 0)
        {
            printf("Os vetores formam uma base para o espaco vetorial.\n");
            adicionarHistorico("Os vetores formam uma base para o espaco vetorial.\n");
        }
        else
        {
            printf("Os vetores NAO formam uma base para o espaco vetorial.\n");
            adicionarHistorico("Os vetores NAO formam uma base para o espaco vetorial.\n");
        }
    }
    else if (dim == 3)
    {
        det = determinante3x3(vetores);
        if (det != 0)
        {
            printf("Os vetores formam uma base para o espaco vetorial.\n");
            adicionarHistorico("Os vetores formam uma base para o espaco vetorial.\n");
        }
        else
        {
            printf("Os vetores NAO formam uma base para o espaco vetorial.\n");
            adicionarHistorico("Os vetores NAO formam uma base para o espaco vetorial.\n");
        }
    }
}

void verificaInjeSobreBije(float a[MAX][MAX + 1], int lin)
{
    int bij = 0;
    int n = lin;
    float temp[MAX][MAX + 1];
    for (int i = 0; i < lin; i++)
    {
        for (int j = 0; j <= lin; j++)
        {
            temp[i][j] = a[i][j];
        }
    }

    eliminacaoDeGauss(temp, n);

    int posto = 0;
    for (int i = 0; i < n; i++)
    {
        int linha_nula = 1;
        for (int j = 0; j < n; j++)
        {
            if (temp[i][j] != 0)
            {
                linha_nula = 0;
                break;
            }
        }
        if (!linha_nula)
        {
            posto++;
        }
    }

    printf("Posto da matriz: %d\n", posto);

    if (posto == n)
    {
        printf("A matriz e Injetiva\n");
        adicionarHistorico("A matriz e Injetiva\n");
        bij = bij + 1;
    }
    else
    {
        printf("A matriz NAO e Injetiva\n");
        adicionarHistorico("A matriz NAO e Injetiva\n");
    }

    if (posto == n)
    {
        printf("A matriz e Sobrejetiva\n");
        adicionarHistorico("A matriz e Sobrejetiva\n");
        bij = bij + 1;
    }
    else
    {
        printf("A matriz NAO e Sobrejetiva\n");
        adicionarHistorico("A matriz NAO e Sobrejetiva\n");
    }

    if (bij == 2)
    {
        printf("A matriz e Bijetiva\n");
        adicionarHistorico("A matriz e Bijetiva\n");
    }
    else
    {
        printf("A matriz NAO e Bijetiva\n");
        adicionarHistorico("A matriz NAO e Bijetiva\n");
    }
}

// Função para multiplicar matrizes 2x2
void multiplicar_matrizes_2x2(float A[2][2], float B[2][2], float res[2][2])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            res[i][j] = 0;
            for (int k = 0; k < 2; k++)
            {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Função para diagonalizar matrizes 2x2
void diagonalizar(float A[2][2])
{
    // Calcular autovalores
    Autovalores autovalores = calcular_autovalores_2x2(A);
    printf("Autovalores:\nlambda1 = %.2f\nlambda2 = %.2f\n", autovalores.lambda1, autovalores.lambda2);

    // Calcular autovetores
    Autovetores autovetores = calcular_autovetores_2x2(A);
    printf("Autovetores:\n");
    printf("Vetor1: [%.2f, %.2f]\n", autovetores.vetor1[0], autovetores.vetor1[1]);
    printf("Vetor2: [%.2f, %.2f]\n", autovetores.vetor2[0], autovetores.vetor2[1]);

    // Matriz de mudança de base P (autovetores)
    float P[2][2] = {
        {autovetores.vetor1[0], autovetores.vetor2[0]},
        {autovetores.vetor1[1], autovetores.vetor2[1]}};

    // Matriz diagonal D (autovalores)
    float D[2][2] = {
        {autovalores.lambda1, 0},
        {0, autovalores.lambda2}};

    // Exibir a matriz diagonalizada correta
    printf("Matriz diagonalizada (D):\n");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%.2f\t", D[i][j]);
        }
        printf("\n");
    }

    printf("Matriz de mudanca de base (P):\n");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%.2f\t", P[i][j]);
        }
        printf("\n");
    }

    // Verificar diagonalização: A = P * D * P^-1
    // Calcular P^-1
    float P_inv[2][2];
    float det = determinante2x2(P);
    if (det == 0)
    {
        printf("A matriz P é singular, não tem inversa.\n");
    }
    P_inv[0][0] = P[1][1] / det;
    P_inv[0][1] = -P[0][1] / det;
    P_inv[1][0] = -P[1][0] / det;
    P_inv[1][1] = P[0][0] / det;

    // Multiplicar P * D
    float PD[2][2];
    multiplicar_matrizes_2x2(P, D, PD);

    // Multiplicar PD * P^-1
    float resultado[2][2];
    multiplicar_matrizes_2x2(PD, P_inv, resultado);

    // Exibir o resultado de A = P * D * P^-1 (deve ser igual a A original)
    printf("Verificacao da diagonalizacao:\n");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%.2f\t", resultado[i][j]);
        }
        printf("\n");
    }

    return 0;
}

// Função para exibir o menu de autovalores e autovetores
void menuAutovaloresAutovetores()
{
    int n;
    printf("Digite a dimensão da matriz (somente 2): ");
    scanf("%d", &n);

    if (n == 2)
    {
        float A[2][2];
        printf("Digite os elementos da matriz 2x2:\n");
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                scanf("%f", &A[i][j]);
            }
        }

        Autovalores autovalores = calcular_autovalores_2x2(A);
        Autovetores autovetores = calcular_autovetores_2x2(A);

        printf("Autovalores:\n");
        printf("Lambda1 = %.2lf\n", autovalores.lambda1);
        printf("Lambda2 = %.2lf\n", autovalores.lambda2);

        printf("Autovetores:\n");
        printf("Vetor1: [%.2f, %.2f]\n", autovetores.vetor1[0], autovetores.vetor1[1]);
        printf("Vetor2: [%.2f, %.2f]\n", autovetores.vetor2[0], autovetores.vetor2[1]);

        adicionarHistorico("Autovalores: Lambda1 = %.2lf, Lambda2 = %.2lf\n", autovalores.lambda1, autovalores.lambda2);
        adicionarHistorico("Autovetores: Vetor1 = [%.2f, %.2f], Vetor2 = [%.2f, %.2f]\n",
                           autovetores.vetor1[0], autovetores.vetor1[1],
                           autovetores.vetor2[0], autovetores.vetor2[1]);
    }
    else
    {
        printf("Dimensão inválida. Somente 2x2 é suportado.\n");
    }
}

// Menu de diagonalização
void menuDiagonalizacao()
{
    int n;
    printf("Digite a dimensão da matriz (somente 2): ");
    scanf("%d", &n);

    if (n != 2)
    {
        printf("Dimensão invalida. Apenas 2x2!\n");
    }

    float A[2][2];
    printf("Digite os elementos da matriz 2x2:\n");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            scanf("%f", &A[i][j]);
        }
    }

    diagonalizar(A);
}

// Função para exibir o menu e obter a escolha do usuário
int exibirMenu()
{
    int escolha;

    printf("=====================================\n");
    printf("    BEM-VINDO AO ALGEBRA CALCULATOR  \n");
    printf("=====================================\n");
    printf("Escolha uma das opcoes abaixo:\n");
    printf("1 - Resolucao de Sistemas Lineares\n");
    printf("2 - Verificacao de Injetividade e Sobrejetividade\n");
    printf("3 - Determinacao de Bases\n");
    printf("4 - Calculo de Autovalores e Autovetores\n");
    printf("5 - Diagonalizar Matriz\n");
    printf("6 - Acessar Historico\n");
    printf("7 - Limpar Historico\n");
    printf("0 - Sair\n");
    printf("=====================================\n");
    printf("Digite sua escolha: ");

    // Lê a escolha do usuário
    scanf("%d", &escolha);

    return escolha;
}

// Função principal
int main()
{
    adicionarHistorico("-------------HISTORICO---------------------\n\n\n");
    int opcao;

    do
    {
        opcao = exibirMenu();

        switch (opcao)
        {
        case 1:
            printf("Opcao 1 escolhida: Resolucao de Sistemas Lineares.\n");
            adicionarHistorico("Opcao 1 escolhida: Resolucao de Sistemas Lineares.\n");
            int n;
            float a[MAX][MAX + 1];
            printf("Digite a dimensao da matriz (2 ou 3): ");
            scanf("%d", &n);
            adicionarHistorico("Dimensao da Matriz: %d\n", n);
            if (n != 2 && n != 3)
            {
                printf("Dimensao invalida. Escolha 2 ou 3.\n");
                adicionarHistorico("Dimensao Invalida\n");
                return 1;
            }

            scanearMatriz(a, n);

            printf("Matriz original:\n");
            imprimirMatriz(a, n);
            adicionarHistorico("Matriz Escalonada:\n", a, n);
            ResultadoSistema resultado = eliminacaoDeGauss(a, n);

            printf("Matriz escalonada:\n");
            imprimirMatriz(a, n);
            adicionarHistorico("Resultado do sistema:\n", a, n);
            if (resultado.sisPossivel)
            {
                resultado = resolucaoSistema(a, n);

                if (resultado.sisDeterminado)
                {
                    printf("Sistema possivel e determinado.\n");
                }
                else
                {
                    printf("Sistema possivel e indeterminado.\n");
                }
            }
            else
            {
                printf("Sistema impossivel.\n");
            }
            break;
        case 2:
            // Chama a função de Verificação de Injetividade e Sobrejetividade
            printf("Opcao 2 escolhida: Verificação de Injetividade, Sobrejetividade e Bijetividade.\n");
            adicionarHistorico("Opcao 2 escolhida: Verificação de Injetividade, Sobrejetividade e Bijetividade.\n");
            printf("Digite a dimensao da matriz (2 ou 3): ");
            scanf("%d", &n);

            if (n != 2 && n != 3)
            {
                adicionarHistorico("Dimensao invalida. Escolha 2 ou 3.\n");
                printf("Dimensao invalida. Escolha 2 ou 3.\n");
                return 1;
            }

            scanearMatriz(a, n);

            printf("Matriz original:\n");
            imprimirMatriz(a, n);
            verificaInjeSobreBije(a, n);
            break;
        case 3:
            // Chama a função de Determinação de Bases
            printf("Opcao 3 escolhida: Determinacao de Bases.\n");
            adicionarHistorico("Opcao 3 escolhida: Determinacao de Bases.\n\n");

            int dim;

            printf("Digite a dimensao dos vetores (2 ou 3): ");
            scanf("%d", &dim);
            adicionarHistorico("Dimensao dos vetores digitadas: %d\n", dim);

            if (dim != 2 && dim != 3)
            {
                adicionarHistorico("Dimensao Invalida\n\n");
                printf("Dimensao invalida. Escolha 2 ou 3.\n");
                return 1;
            }

            float vetores[MAX][MAX];

            printf("Digite os vetores:\n");
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    scanf("%f", &vetores[i][j]);
                }
            }

            verificarBase(vetores, dim);
            break;
        case 4:
            printf("Opcao 4 escolhida: Cálculo de Autovalores e Autovetores.\n");
            menuAutovaloresAutovetores();
            adicionarHistorico("Opcao 4 escolhida: Cálculo de Autovalores e Autovetores.\n");
            break;
        case 5:
            printf("Opcao 5 escolhida: Diagonalizar matriz.\n");
            menuDiagonalizacao();
            adicionarHistorico("Opcao 5 escolhida: Diagonalizar matriz.\n");
            break;
        case 6:
            imprimirHistorico();
            printf("---------------FIM HISTORICO---------------------\n\n\n");
            break;
        case 7:
            limparHistorico();
            adicionarHistorico("-------------------INICIO HISTORICO------------------\n\n\n\n");
            printf("---------------LIMPAR HISTORICO---------------------\n\n\n");
            break;
        case 0:
            adicionarHistorico("Programa Finalizado\n\n\n\n");
            printf("Saindo do programa...\n");
            break;
        default:
            printf("Opcao invalida! Por favor, escolha uma opção valida.\n");
        }

        printf("\n");

    } while (opcao != 0);

    return 0;
}
