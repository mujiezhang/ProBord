#!/bin/bash
# Usage: bash ./prepare_blast_db.sh [options] <genus_name> [group]

show_help() {
  echo "Usage: bash $0 [options] <genus_name> [bacteria|archaea]"
  echo "Downloads NCBI genomes and builds BLAST database for specified genus"
  echo ""
  echo "Options:"
  echo "  -h, --help      Show this help message"
  echo "  -t, --threads N Set number of download threads (default: 4)"
  echo "  -o, --output DIR Specify custom output directory (default: use genus name)"
  echo ""
  echo "Examples:"
  echo "bash $0 Mannheimia bacteria"
  echo "bash $0 Methanobrevibacter archaea"
  exit 0
}

# 默认参数
THREADS=4
GROUP="bacteria"
OUTPUT_DIR=""

# 解析参数
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      show_help
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    *)
      # 第一个非选项参数为属名
      if [[ -z "$GENUS" ]]; then
        GENUS="$1"
      # 第二个非选项参数为分组
      elif [[ -z "$GROUP_ARG" ]]; then
        GROUP_ARG="$1"
      else
        echo "Error: Too many arguments"
        show_help
        exit 1
      fi
      shift
      ;;
  esac
done

# 验证必填参数
if [[ -z "$GENUS" ]]; then
  echo "Error: Genus name is required"
  show_help
  exit 1
fi

# 使用自定义分组（如果提供）
if [[ -n "$GROUP_ARG" ]]; then
  GROUP="$GROUP_ARG"
fi

# 设置输出目录 (修改点：默认直接使用属名)
if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="$GENUS"  # 直接使用属名作为目录名
fi

# 创建数据库目录
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" || exit 1

echo "======================================================"
echo "Building BLAST database for: $GENUS ($GROUP)"
echo "Output directory: $PWD"
echo "Download threads: $THREADS"
echo "======================================================"

# 下载基因组
echo "[1/4] Downloading $GENUS genomes from $GROUP..."
ncbi-genome-download -F fasta --genera "$GENUS" -o "$GENUS" --flat-output -P -p "$THREADS" "$GROUP"

# 检查下载结果
if [[ ! "$(ls -A "$GENUS")" ]]; then
  echo "Error: No genomes downloaded. Possible reasons:" >&2
  echo "  - Genus name misspelled (case-sensitive)" >&2
  echo "  - No genomes available for $GENUS in $GROUP" >&2
  echo "  - NCBI API connection issue" >&2
  exit 1
fi

# 处理基因组
echo "[2/4] Processing downloaded genomes..."
gunzip "$GENUS"/*.gz 2>/dev/null
cat "$GENUS"/*.fna > "${GENUS}.fna"

# 创建BLAST数据库
echo "[3/4] Building BLAST database (this may take several minutes)..."
makeblastdb -in "${GENUS}.fna" -dbtype nucl -out "$GENUS" -title "${GENUS}_DB"

# 清理临时文件
echo "[4/4] Cleaning temporary files..."
rm -rf "$GENUS"
rm "${GENUS}.fna"

# 完成信息
echo "======================================================"
echo "Successfully created BLAST database for $GENUS"
