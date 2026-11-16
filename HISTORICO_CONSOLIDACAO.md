# Consolidacao do sistema CPC — 2026-11-16

Registro do que foi feito ao unificar as copias que existiam em paralelo.

## Situacao anterior

| Local | Conteudo | Git |
|---|---|---|
| `<raiz-do-repositorio>` | 64 arquivos, 14 MB — CLIMA/, MAPAS_CPC/, NFE/, SAE/, SIN/ + scripts de raiz | nao |
| `<pasta-antiga>` | 14 arquivos, 338 KB — OPERACIONAL/, PLOTAGEM/, SAE/, SIN/ | sim, `github.com/deepblueisdark/CPC` |
| `<servidor>` (conta `simop`) | versao operacional | nao verificada |

## Problemas encontrados

1. **Pipeline partido entre as duas copias.** `roda_SAE_opera.sh`
   (so no Devel) chamava `cpc_gera_plota_chuva_acumulada.py` (so no
   GitHub) e `cpc_plota_clima_unificado.py` (so no Devel). Nenhuma das
   duas copias rodava sozinha.
2. **Dois `cpc.config` incompativeis.** O do GitHub tinha os caminhos
   entre aspas literais (`SHAPEFILE_DIR = "<caminho-wsl>/..."`), o que quebrava
   os `os.path.join`, e usava caminhos WSL. O do Devel estava correto
   para Windows, mas tinha perdido `MAPA_DIR` e `saida_excel`.
3. **`roda_SAE_opera.sh` nao rodava com nenhum dos dois configs:**
   pedia `--regiao SAE` (so existiam AMS, BRASIL e SAARA) e
   `--cores escala_cores_mensal` (chave inexistente). Os scripts
   levantam `ValueError` nos dois casos.
4. **`.gitignore` ignorava `cpc.config`**, que estava rastreado no
   repositorio — regra inocua e enganosa.
5. **Devel nao era repositorio git**, apesar de ser onde o
   desenvolvimento acontecia.
6. **Tres versoes de `cpc_plota_clima`** divergentes (V1 no Devel, V1
   diferente no GitHub, V2), alem do `unificado` que o proprio
   cabecalho declara substituir V1 e V2.
7. **Codigo misturado com saida:** 42 PNGs em `MAPAS_CPC/`, 4 PNGs em
   `CLIMA/`, planilha de 7 MB em `NFE/`.
8. **Todos os scripts procuravam `cpc.config` no diretorio atual**, o
   que impedia qualquer organizacao em subpastas.
9. **Referencias ausentes:** `SAE_TODAS_AS_BACIAS_HICON_POLY.shp` e
   quatro planilhas usadas pelas rotinas SIMOP.

Os arquivos em comum (`cpc_download.py`, `cpc_calcbacia.py`,
`SAE/cpc_SIMOP_SAE.py`, `SIN/cpc_SIMOP_SIN.py`,
`SIN/cpc_SIMOP_SIN_GRAFICOS.py`) estavam **byte a byte identicos** nas
duas copias — nao houve conflito de conteudo.

## O que foi feito

### Backups (antes de qualquer mudanca)
- `_backup/BACKUP_Devel_CPC_2026-11-16.tar.gz` (13 MB)
- `_backup/BACKUP_GITHUB_CPC_2026-11-16.tar.gz` (113 KB, inclui o `.git` completo)

### Scripts trazidos do GitHub (existiam so la)
- `PLOTAGEM/cpc_gera_plota_chuva_acumulada.py`
- `PLOTAGEM/cpc_plota.py`
- `PLOTAGEM/cpc_climatology.py`

### Reorganizacao
- `cpc_download.py`, `cpc_calcbacia.py` -> `OPERACIONAL/`
- `cpc_plota_clima_unificado.py` -> `PLOTAGEM/`
- `roda_SAE_opera.sh` -> `scripts/`
- 49 PNGs -> `SAIDAS/MAPAS_CPC/` e `SAIDAS/CLIMA/`
- `chuva_diaria_CPC.xlsx`, `chuva_NFE.xlsx` -> `SAIDAS/PLANILHAS/`
- `chuva_diaria_CPC_NFE.xlsx` (entrada, 7 MB) -> `DADOS/`
- pastas `CLIMA/` e `MAPAS_CPC/` removidas (vazias)

### Arquivadas em `_legado/`
- `cpc_plota_clima_V1_devel.py` (era `CLIMA/cpc_plota_clima.py`)
- `cpc_plota_clima_V1_github.py` (era `PLOTAGEM/cpc_plota_clima.py`)
- `cpc_plota_clima_V2.py`

Ativo: apenas `PLOTAGEM/cpc_plota_clima_unificado.py`.

### Correcoes de codigo
- **Resolucao de caminhos** em 6 scripts (`cpc_download.py`,
  `cpc_calcbacia.py`, `cpc_plota_clima_unificado.py`,
  `cpc_gera_plota_chuva_acumulada.py`, `cpc_plota.py`,
  `cpc_climatology.py`): adicionado `_raiz_projeto()`, que procura o
  `cpc.config` na pasta do script e acima. Agora rodam de qualquer
  diretorio.
- `MAPA_DIR` em `cpc_calcbacia.py` passou de `./MAPAS_CPC` para
  `SAIDAS/MAPAS_CPC` na raiz do projeto.
- `roda_SAE_opera.sh` reescrito: resolve a raiz via `BASH_SOURCE`,
  aceita a variavel `SHAPEFILE`, falha com mensagem clara se o
  shapefile nao existir, e grava em `SAIDAS/CLIMA/`.

### `cpc.config` unificado
- Caminhos Windows ativos; equivalentes Linux/WSL comentados.
- Aspas literais dos caminhos removidas.
- `mapa_dir` e `saida_excel` recuperados, apontando para `SAIDAS/`.
- Adicionada a chave `escala_cores_mensal`.
- Adicionada a secao `[SAE]`.

### Git
- Repositorio novo iniciado nesta pasta (historico anterior preservado
  apenas no `.tar.gz` de backup e no GitHub remoto).
- `.gitignore` reescrito: ignora `SAIDAS/`, `DADOS/`, `_backup/`,
  `__pycache__/`, `*.nc` e `desktop.ini`; a regra inocua sobre
  `cpc.config` foi removida e o arquivo segue rastreado.

### Pasta antiga
- Conteudo de `<pasta-antiga>` movido para
  `<pasta-antiga>\_to_delete\`, para conferencia e exclusao
  manual.

## Nao feito

- **Comparacao com o servidor** `<servidor>` — esta sessao nao
  alcanca o host (SSH fora da allowlist de rede).
- **Push para o GitHub** — o historico novo divergiu do remoto; um push
  seria rejeitado. Decida entre `--force` (sobrescreve o remoto) ou um
  branch novo.
- Os valores de `[SAE]` e `escala_cores_mensal` precisam da sua
  confirmacao.

---

# Segunda rodada — comparacao com o servidor (2026-11-16, tarde)

A sessao nao alcanca `<servidor>` por rede (nenhuma porta sai da
allowlist) e o computer use recusa digitacao em IDEs e terminais
(`tier: "click"`). A coleta foi feita pelo WSL da maquina local, com
`scripts/coleta_servidor.sh`, que alcanca o servidor e enxerga o `E:`.

## O que havia no servidor

`/home/<usuario>/PRODUCAO/CPC` tem **tres** copias internas:

| Diretorio | Papel |
|---|---|
| `CPC/` (raiz) | copia antiga, rastreada no git, 8 commits atras do GitHub |
| `CPC/OPERA/` | **producao** — e o que o cron executa; nao rastreada |
| `CPC/SAE/` | sandbox de desenvolvimento; nao rastreada |

O git da raiz estava em `624e45d`, oito commits atras do
`origin/main`, com `cpc.config` modificado e `OPERA/`, `SAE/`,
`cpc.config.local` e `log_cpc.txt` nao rastreados. Ou seja: o que
roda em producao nunca esteve sob versionamento.

## Resultado da comparacao

Comparacao por hash, normalizando CRLF e espacos finais.

**Todo o codigo que o servidor executa em `OPERA/` e identico ao que
foi consolidado** — os 10 scripts batem byte a byte. A consolidacao
da manha acertou o codigo que importa.

Faltavam duas pecas, ambas so no servidor:

1. **`OPERA/rodar_cpc.sh`** — o orquestrador do cron, que amarra o
   sistema inteiro. Nao existia em nenhuma copia local. Portado para
   `scripts/rodar_cpc.sh`; original em `_legado/servidor/`.
2. **`OPERA/cpc.config`** — o config de producao, com as escalas e a
   regiao `[SAE]` reais.

E tres scripts de funcionalidade ausente, da sandbox `SAE/`:

- `cpc_calcula_anomalia.py` -> `PLOTAGEM/` (introduz dependencia `xarray`)
- `cpc_plota_anomalia.py` -> `PLOTAGEM/`
- `cpc_plota_chuva_acumulada.py` e `cpc_calcula_anomalia_old.py` -> `_legado/`

## Correcoes que a comparacao obrigou

- **Os valores que eu havia chutado estavam errados.** A
  `escala_cores_mensal` de producao tem 31 niveis de 0.1 a 300 mm/mes;
  eu havia posto 15 niveis de 1 a 1000. A regiao `[SAE]` real e
  lat -8/-21, lon 285/305; eu havia posto -21/-6 e 284/306.
  Ambos substituidos pelos valores de producao.
- **Escalas que faltavam** foram acrescentadas: `escala_diaria`,
  `escala_clima`, `escala_clima_simples`, e as secoes `[ANOMALIA]` e
  `[ANOMALIA_PER]` (estas vindas da sandbox `SAE/cpc.config`).
- **`[SAE]` tem `min_lat` maior que `max_lat`** (-8 contra -21),
  invertido em relacao as outras secoes. Copiado verbatim porque e
  assim que roda em producao; sinalizado no `cpc.config` e no README.
- **As duas "versoes divergentes" de `cpc_plota_clima.py` eram o mesmo
  arquivo.** Diferiam apenas em CRLF. Os dois arquivos que eu havia
  arquivado como `_V1_devel` e `_V1_github` foram colapsados em
  `_legado/cpc_plota_clima_V1.py`.
- **`DADOS/` foi substituido por `OPERA/`.** A leitura dos scripts
  SIMOP mostrou que eles leem e gravam as planilhas acumuladoras por
  nome simples no diretorio corrente, entao o diretorio de trabalho
  precisa ser plano, como no servidor. As tres planilhas locais foram
  movidas para `OPERA/`.
- **O shapefile deixou de ser pendencia**: `SAE_TODAS_AS_BACIAS_HICON_POLY`
  veio do servidor para `SHAPES/SAE/`, e o `roda_SAE_opera.sh` passou a
  apontar para ele por padrao.
- `.gitattributes` ganhou `*.sh text eol=lf`, ja que os shells rodam
  no servidor Linux.

## O que o servidor tem e nao foi trazido

- As quatro planilhas acumuladoras (12 MB), excluidas do pacote de
  coleta: `CHUVA_OBSERVADA.xlsx`, `CHUVA_OBSERVADA_CPC_COMPLETA.xlsx`,
  `CHUVA_CPC_COMPLETA_SIN.xlsx`, `chuva_SIN.xlsx`. Sem elas a cadeia
  SIMOP nao roda localmente.
- Os shapes de SIN e NFE (`/home/<usuario>/STUFF/SHAPES/{SIN,NFE}`).
- As copias divergentes e mais antigas da sandbox `SAE/`
  (`cpc_calcbacia.py`, `roda_SAE_opera.sh`, `cpc.config`) — eram
  rascunhos nunca promovidos a producao. Os dois `cpc.config` do
  servidor estao em `_legado/servidor/` para referencia.

## Ainda em aberto

- O servidor continua rodando a partir de `OPERA/`, fora do git. Nada
  foi alterado la. Decidir se producao passa a apontar para o
  repositorio consolidado, e como.
- Push para o GitHub: o historico novo divergiu do remoto.

---

# Terceira rodada — repositorio de codigo (2026-11-16, fim da tarde)

Decisao do usuario: o repositorio guarda **todos os scripts** das tres
aplicacoes (SAE, SIN, NFE) e **nenhum operacional**. O que roda no dia
a dia e a producao no servidor; a pasta local e so desenvolvimento.

## Acervo de scripts completado

Acrescentados a `_legado/servidor/` os quatro arquivos do servidor que
ainda nao tinham representacao: `cpc_calcbacia.py.sandbox_SAE`,
`cpc_calcbacia.py.raiz`, `roda_SAE_opera.sh.sandbox_SAE` e
`cpc.config.local`.

Verificacao por hash, desconsiderando CRLF e o bloco de resolucao de
caminhos: **os 33 scripts e configs do servidor estao todos no
repositorio**, como versao ativa ou copia verbatim. As unicas
diferencas em relacao ao servidor sao os patches intencionais da
consolidacao.

## Operacional retirado da pasta local

- 3 planilhas de `OPERA/` (7,2 MB) e 49 PNGs de `SAIDAS/` (6,5 MB)
  apagados. As duas pastas ficam vazias, so como ponto de montagem da
  execucao local.
- `SHAPES/` removido do repositorio, inclusive o shapefile SAE que
  havia sido trazido na segunda rodada. Shapefile e dado operacional.
  O `roda_SAE_opera.sh` voltou a apontar para o caminho do servidor,
  sobreponivel pela variavel `SHAPEFILE`.
- `.gitignore` passou a barrar `*.xlsx`, `SHAPES/` e as extensoes de
  shapefile, para nao entrarem por descuido.

## Backups movidos para fora do repositorio

De `CPC/_backup/` para `<pasta-de-backups>\`,
conferidos por md5 e com os tar.gz testados antes de apagar os
originais:

- `BACKUP_Devel_CPC_2026-11-16.tar.gz` (13 MB)
- `BACKUP_GITHUB_CPC_2026-11-16.tar.gz` (113 KB, com o .git original)
- `CPC_SERVIDOR.tar.gz` (796 KB)
- `INVENTARIO_SERVIDOR.txt` (30 KB)

Tudo que foi apagado hoje esta nesses arquivos.

## Pendencia tecnica

`__pycache__/cpc_calcbacia.cpython-312.pyc` nao pode ser removido:
`Permission denied`, provavelmente trava do OneDrive. Esta ignorado
pelo git; apagar na mao quando o OneDrive liberar.

---

# Quarta rodada — sincronizacao com a producao (2026-11-16, noite)

Decisao do usuario: a producao **nao** muda. O servidor continua
rodando de `~/PRODUCAO/CPC/OPERA/`, fora do git. O que o repositorio
passa a fazer e **capturar** esse diretorio, para o codigo que roda de
fato ficar sob historico de versao.

## Sincronizacao de captura, em um sentido so

Criado `scripts/sincroniza_opera.sh`, que roda no WSL da maquina de
desenvolvimento — a unica ponte que alcanca o servidor por SSH e ao
mesmo tempo enxerga o `E:`. O script:

1. le a arvore completa do `OPERA/` remoto para `_producao/ARVORE_OPERA.txt`,
   com tamanhos, datas e hashes, inclusive dos arquivos que nao sao
   copiados;
2. traz codigo e configs para `_producao/OPERA/`, substituindo o
   snapshot inteiro, para que arquivo apagado no servidor desapareca
   tambem daqui;
3. gera `_producao/MANIFESTO.txt` com o md5 de cada arquivo;
4. mostra **o que mudou na producao desde o ultimo snapshot commitado**
   — a deteccao de drift que interessa no dia a dia;
5. compara o snapshot com o codigo ativo do repositorio, normalizando
   CRLF e espaco no fim da linha.

Acesso remoto so de leitura (`ssh` + `tar` para stdout). O script nunca
escreve no servidor.

Dado operacional fica fora do snapshot, como no resto do projeto.
`log_cpc.txt` precisou de exclusao por nome: casa com `*.txt` e teria
entrado junto com o codigo.

## `_producao/`

- `README.md` — o que a pasta e, como atualizar, o que entra e o que nao.
- `DIVERGENCIAS.md` — as diferencas **intencionais** entre repositorio e
  producao, para a etapa 5 ter como ser lida. Sem essa lista, sete
  arquivos aparecem como divergentes e nao da para distinguir patch de
  drift.

Na conferencia, o levantamento da terceira rodada estava desatualizado:
sao **sete** os scripts com o bloco `_raiz_projeto()`, nao seis —
`PLOTAGEM/cpc_plota_anomalia.py` tambem o tem. `cpc_calcula_anomalia.py`
nao entra na lista porque nao le o `cpc.config`: recebe tudo por
`argparse`.

## `[SAE]` corrigido no repositorio, nao na producao

As latitudes da secao `[SAE]` passaram a `min_lat = -21` e
`max_lat = -8`, na mesma convencao das outras secoes. A producao segue
com `-8` / `-21`, como esta desde 2025-10.

E uma divergencia deliberada, registrada no `cpc.config`, no
`README.md` e em `_producao/DIVERGENCIAS.md`. **Antes de levar o config
ao servidor, gerar uma imagem do SAE com os valores novos e comparar
com a do dia** — o recorte do mapa pode mudar, e por enquanto o que
esta no ar e a versao invertida.

## Nao feito

- **Push para o GitHub**, por decisao do usuario nesta rodada. O
  historico local segue divergente de `origin/main` (`ec25759`).
- **Primeira execucao do `sincroniza_opera.sh`**: depende de rodar no
  WSL. Enquanto nao rodar, `_producao/OPERA/` esta vazia.

---

# Quinta rodada — rotinas de empresa fora do repositorio (2026-11-16, fim do dia)

Correcao de rumo. O escopo registrado na terceira rodada — "o
repositorio guarda **todos os scripts** das tres aplicacoes (SAE, SIN,
NFE)" — estava errado.

**Decisao dele:** o repositorio guarda **apenas os scripts
principais**. As rotinas especificas de cada empresa ficam **so no
servidor**, porque sao proprias de cada uma. (O caso ABCLUP e manual,
nao tem script.)

## O que saiu

As cinco rotinas SIMOP, que estavam em **dois** lugares:

| Onde | Desde |
|---|---|
| `SAE/`, `SIN/`, `NFE/` | `a019e0e`, o **primeiro** commit |
| `_producao/OPERA/` | `84f8bda`, o snapshot da quarta rodada |

Por estarem no primeiro commit, estavam em todos os dez. `git rm` nao
bastava: foram purgadas do historico, o que so foi limpo porque nada
havia sido pushado.

O snapshot caiu de 12 para 7 arquivos. Consequencia que vale registrar:
**os cinco arquivos que apareciam como identicos ao codigo ativo eram
justamente as rotinas de empresa.** Dos 7 que restam, nenhum e
identico — 4 diferem pelo patch de caminhos, 2 por serem versao
portada e 1 e o `cpc.config`.

## O que impede a volta

O filtro mora no proprio `scripts/sincroniza_opera.sh`, na variavel
`EXCLUI_EMPRESA` (`*SIMOP*`), e vale para as **tres** leituras do
servidor: a arvore, a lista de hashes e a captura. Assim nem o conteudo
nem os nomes voltam. A etapa 1 passou a reportar a **quantidade**
encontrada, para que uma rotina nova no servidor nao passe em branco.

Testado contra um diretorio que imita o `OPERA/`: as cinco rotinas e o
`log_cpc.txt` ficam fora, os scripts principais entram.

## `rodar_cpc.sh` com as rotinas fora

O orquestrador chamava as rotinas por `$RAIZ/SAE/...`, que nao existe
mais. Agora usa a variavel `EMPRESAS`, que por padrao aponta para o
diretorio de trabalho — onde elas moram em producao — e confere as
cinco **antes** de comecar, parando com mensagem clara se faltar
alguma. Conferido: falha com codigo 1 sem criar lock, log nem arquivo.

## O que continua visivel, e nao foi escondido

Saiu o **codigo** das rotinas. Continuam no repositorio os **nomes** e
a **sequencia de chamadas**, porque aparecem no orquestrador, no
`cpc.config`, nos `_legado/servidor/` e na descricao do pipeline no
README. Tirar isso tambem seria outra decisao, nao foi pedida e nao foi
feita.

## Registros limpos

`_producao/MANIFESTO.txt` perdeu as 5 linhas de hash e
`_producao/ARVORE_OPERA.txt` as 10 entradas correspondentes, com nota
dizendo que foram capturados antes do filtro existir. As capturas
seguintes ja saem filtradas na origem.

## Backup antes da purga

`_backup/ANTES_DA_PURGA_2026-11-16_1805.tar.gz` (345 KB), com o `.git`
nos dez commits originais. Testado: restaura repositorio funcional,
`fsck` limpo, as 10 entradas SIMOP preservadas. **Esse tarball e o
unico lugar onde o historico antigo existe** — e esta dentro de
`_backup/`, que o `.gitignore` barra. Mova para fora da pasta do
repositorio.

---

# Sexta rodada — preparar o repositorio para ser publico (2026-11-16, noite)

Ao decidir que o repositorio no GitHub passaria a espelhar esta pasta,
veio a descoberta que mudou o plano: **o repositorio publico nao expunha
nada do servidor, e o trabalho do dia tinha acrescentado tudo.**

O publico (`ec25759`, de 2025-11-05) tinha zero ocorrencias de IP, de
usuario SSH, de hostname ou de caminho `/home/...`; so caminhos WSL da
maquina local. Esta pasta, depois do dia, tinha o IP do servidor, o
usuario, o hostname, o layout `/home/<usuario>/{PRODUCAO,STUFF,DADOS}` e
o horario do cron em **20 arquivos rastreados** — nos scripts de coleta
e sincronizacao e na documentacao.

Publicar assim diria onde a producao esta, com que usuario entrar, que
diretorios existem e quando os jobs rodam.

## Decisao: publico, mas sem identidade de ambiente

Nada que diga onde a producao roda fica rastreado. O que saiu do git:

| Saiu | Por que |
|---|---|
| `_producao/` | espelho verbatim, carrega a config real da producao |
| `_legado/servidor/` | copias verbatim do servidor, com caminhos reais |
| `cpc.config` | caminhos desta maquina nas chaves ativas |

Os tres continuam **em disco**, so deixaram de ser versionados.

## O que entrou no lugar

Dois modelos rastreados, e o real fica fora:

- `cpc.local.example` -> `cpc.local`: servidor, `OPERA_REMOTO` e os
  caminhos dos shapes. Usa a forma `${VAR:=valor}`, que so define se a
  variavel ainda nao existir — assim um valor passado pelo ambiente
  continua tendo precedencia.
- `cpc.config.example` -> `cpc.config`: as escalas e as regioes (que
  sao o conteudo que vale compartilhar) com os dois caminhos de
  `[DEFAULT]` como placeholder.

Os quatro scripts de shell passaram a carregar o `cpc.local` e a parar
com mensagem clara se faltar variavel. Nenhum tem mais valor real
embutido.

**Detalhe que vale saber:** o `cpc.config` nao e opcional num clone
novo. O bloco `_raiz_projeto()` localiza a raiz do projeto procurando
justamente esse arquivo; sem ele, cada script cai no diretorio
corrente. Por isso "Primeiros passos" no README manda copiar os dois
modelos antes de qualquer coisa.

## Deteccao de drift sem git

A etapa 4 do `sincroniza_opera.sh` usava `git status` em
`_producao/OPERA`, que deixou de ser rastreado. Agora compara o
manifesto novo com o anterior, guardado ao lado em
`MANIFESTO_ANTERIOR.txt`, e classifica cada arquivo como ALTERADO,
NOVO ou SUMIU. Testado nos tres cenarios e no caso do primeiro
snapshot, sem anterior com que comparar.

A perda real: o historico do que rodou em producao deixa de existir no
git. Fica o manifesto anterior, que cobre a comparacao de uma
sincronizacao para a seguinte, mas nao a serie inteira.

## Limpeza dos arquivos que ficaram

IP, usuario, hostname, caminhos do servidor, caminhos locais da maquina
e o horario do cron substituidos por placeholders em `README.md`,
`HISTORICO_CONSOLIDACAO.md`, `cpc.config.example` e nos quatro scripts.
Cuidado tomado: as substituicoes casam `<usuario>@` e `/home/<usuario>`, nunca
`simop` solto, para nao estragar as mencoes a `SIMOP` em maiusculas,
que sao outra coisa.

Verificado: zero ocorrencias de IP, usuario SSH, hostname, caminhos do
servidor ou caminhos de maquina em qualquer arquivo rastreado. As
mencoes a "OneDrive" que sobraram sao sobre a trava de arquivo do
OneDrive, nao caminhos.

## GitHub

O repositorio publico vai ser **apagado e recriado**, para que nenhum
objeto antigo fique acessivel por SHA — um `force-push` trocaria o topo
mas deixaria os commits antigos recuperaveis. Nele estavam expostas,
desde 2025-11, tres rotinas de empresa: `cpc_SIMOP_SAE.py`,
`cpc_SIMOP_SIN.py` e `cpc_SIMOP_SIN_GRAFICOS.py`.

O que ja foi clonado ou forkado nesse periodo nao tem volta.
