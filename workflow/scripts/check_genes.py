import sys
import os
import re
from Bio import SeqIO

def parse_coverage(header):
    """Tenta extrair a cobertura do cabeçalho (especial para SPAdes/GetOrganelle)."""
    match = re.search(r"cov_([\d\.]+)", header)
    if match:
        return float(match.group(1))
    return None

def parse_length_trinity(header):
    """Tenta extrair o tamanho do cabeçalho do Trinity (len=X)."""
    match = re.search(r"len=(\d+)", header)
    if match:
        return int(match.group(1))
    return None

def check_genes_smart(input_fasta, report_file, ref_fasta=None):
    print(f"A processar o ficheiro: {input_fasta}")
    
    # -------------------------------------------------------------------------
    # MODO B: Se o utilizador passou uma referência válida via linha de comandos
    # -------------------------------------------------------------------------
    if ref_fasta and os.path.exists(ref_fasta) and os.path.getsize(ref_fasta) > 0:
        print(f"-> [MODO B] Referência detetada ({ref_fasta}). A executar validação por BLAST...")
        
        blast_cmd = f"tblastn -query {ref_fasta} -subject {input_fasta} -outfmt '6 qseqid sseqid pident length evalue' -evalue 1e-5"
        try:
            import subprocess
            result = subprocess.run(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
            blast_output = result.stdout.strip()
            
            with open(report_file, "w") as f:
                f.write(f"MitoCatch - Relatório via Referência Proteica (BLAST)\n Analisado: {os.path.basename(input_fasta)}\n\n")
                if blast_output:
                    f.write(f"{'GENE_REF':<15} {'CONTIG':<30} {'%_IDENTIDADE':<15} {'ALINHAMENTO':<15}\n")
                    f.write("-" * 75 + "\n")
                    for line in blast_output.split("\n"):
                        p = line.split("\t")
                        f.write(f"{p[0]:<15} {p[1]:<30} {p[2]:<15} {p[3]:<15}\n")
                else:
                    f.write("Aviso: Nenhum gene da referência foi encontrado nesta montagem.\n")
            print(f"Sucesso: Relatório de BLAST gravado em {report_file}")
            return
        except Exception as e:
            print(f"Erro ao executar BLAST, a reverter para Modo Automático (Matemático)... Erro: {e}")

    # -------------------------------------------------------------------------
    # MODO C (Por omissão): Filtragem Matemática por Tamanho e Cobertura
    # -------------------------------------------------------------------------
    print("-> [MODO C] Sem referência. A filtrar contigs por tamanho e cobertura...")
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    filtered_records = []
    
    for record in records:
        header = record.description
        seq_len = len(record.seq)
        
        # Se for Trinity, tenta apanhar o tamanho pelo cabeçalho, senão usa o real
        trinity_len = parse_length_trinity(header)
        actual_len = trinity_len if trinity_len else seq_len
        
        cov = parse_coverage(header)
        
        # Regra de Filtragem:
        # 1. Se tiver cobertura (GetOrganelle), guarda se cov >= 10.0 E tamanho >= 500bp
        # 2. Se não tiver cobertura (Trinity), guarda apenas se o tamanho >= 500bp
        if cov is not None:
            if cov >= 10.0 and actual_len >= 500:
                filtered_records.append((record.id, actual_len, f"{cov:.2f}x"))
        else:
            if actual_len >= 500:
                filtered_records.append((record.id, actual_len, "N/A (Trinity)"))

    # Escrever o relatório final limpo
    with open(report_file, "w") as f:
        f.write("=====================================================================\n")
        f.write("MitoCatch - Relatório de Filtragem Automática (Filtro de Ruído)\n")
        f.write(f"Análise do ficheiro: {os.path.basename(input_fasta)}\n")
        f.write("=====================================================================\n\n")
        f.write(f"Total de contigs analisados na montagem bruta: {len(records)}\n")
        f.write(f"Total de contigs guardados (filtros aplicados): {len(filtered_records)}\n\n")
        
        if filtered_records:
            f.write(f"{'ID_CONTIG':<40} {'TAMANHO (bp)':<15} {'COBERTURA':<15}\n")
            f.write("-" * 75 + "\n")
            # Ordenar os contigs do maior para o menor para ficar arrumado
            for cid, clen, ccov in sorted(filtered_records, key=lambda x: x[1], reverse=True):
                f.write(f"{cid:<40} {clen:<15} {ccov:<15}\n")
        else:
            f.write("Aviso: Nenhum contig passou os critérios mínimos (Tamanho >= 500bp, Cobertura >= 10x).\n")

    print(f"Sucesso: Relatório automático guardado em: {report_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python check_genes.py <input_fasta> <output_report_txt> [opcional_ref_fasta]")
        sys.exit(1)
        
    # Se o Snakemake passar um terceiro argumento, será a referência
    optional_ref = sys.argv[3] if len(sys.argv) > 3 else None
    check_genes_smart(sys.argv[1], sys.argv[2], optional_ref)
