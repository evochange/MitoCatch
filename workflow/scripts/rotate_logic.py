import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

def rotate_with_pure_python(input_file, output_file):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ref_protein_file = os.path.join(script_dir, "cox1_ref.protein.fasta")

    if not os.path.exists(ref_protein_file):
        print(f"Erro: Não encontrei a referência proteica em {ref_protein_file}")
        sys.exit(1)

    # 1. Carregar a proteína de referência cox1
    ref_record = next(SeqIO.parse(ref_protein_file, "fasta"))
    ref_peptide = str(ref_record.seq).upper()
    # Usar um "pedaço" conservado da proteína (ex: os primeiros 15 aminoácidos) para garantir match
    query_subseq = ref_peptide[:15] 

    # 2. Ler TODOS os contigs da montagem
    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        print(f"Erro: Ficheiro {input_file} vazio.")
        sys.exit(1)

    print(f"Carregados {len(records)} contigs do ficheiro de entrada.")

    id_contig_rodado = None
    posicao_rodagem = 0

    # 3. Procurar o gene traduzindo os contigs nos 6 frames (3 forward + 3 reverse complement)
    print("A procurar o gene COX1 nos contigs...")
    for record in records:
        nuc_seq = record.seq
        nuc_len = len(nuc_seq)
        
        # Só analisar contigs com tamanho suficiente para albergar o gene
        if nuc_len < 100:
            continue

        found = False
        # Testar as duas strands: Forward (1) e Reverse (-1)
        for strand, seq_atua in [(1, nuc_seq), (-1, nuc_seq.reverse_complement())]:
            # Testar os 3 frames de leitura de codões
            for frame in range(3):
                try:
                    # Traduzir para proteína (ignora o aviso de tamanho não múltiplo de 3)
                    translated_protein = str(seq_atua[frame:].translate(to_stop=False)).upper()
                    
                    if query_subseq in translated_protein:
                        amino_pos = translated_protein.index(query_subseq)
                        # Converter a posição de aminoácidos de volta para nucleótidos
                        nuc_pos_in_strand = (amino_pos * 3) + frame
                        
                        if strand == 1:
                            posicao_rodagem = nuc_pos_in_strand
                        else:
                            # Se estiver na strand reverse, calcula a posição real na strand forward
                            posicao_rodagem = nuc_len - nuc_pos_in_strand
                        
                        id_contig_rodado = record.id
                        found = True
                        break
                except Exception:
                    continue
            if found:
                break
        if found:
            break

    # 4. Se encontrou, rodar o contig alvo
    if id_contig_rodado:
        print(f"Gene COX1 detetado no contig [{id_contig_rodado}] próximo da posição {posicao_rodagem}. A rodar contig...")
        for record in records:
            if record.id == id_contig_rodado:
                full_seq = str(record.seq)
                rotated_seq = full_seq[posicao_rodagem:] + full_seq[:posicao_rodagem]
                record.seq = Seq(rotated_seq)
                print(f"Contig [{record.id}] rodado com sucesso!")
                break
    else:
        print("Aviso: Não foi possível mapear o COX1 via tradução proteica. Gravado sem rotação.")

    # 5. Gravar TODOS os contigs no ficheiro de output
    SeqIO.write(records, output_file, "fasta")
    print(f"Ficheiro completo guardado com sucesso em: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python rotate_logic.py <input_fasta> <output_fasta>")
        sys.exit(1)

    rotate_with_pure_python(sys.argv[1], sys.argv[2])

