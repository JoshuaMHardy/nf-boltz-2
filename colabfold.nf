// main.nf

// Workflow

workflow {
    // Prepare sequences
    PREPARE_SEQUENCES()
    
    // Only proceed if there are protein sequences needing MSA
    PREPARE_SEQUENCES.out.seq_info
        .filter { pred_name, seq_json, chains -> 
            chains.text.trim().length() > 0 
        }
        .set { sequences_to_process }
    
    // Submit to ColabFold (unpaired + optionally paired)
    SUBMIT_COLABFOLD(sequences_to_process)
    
    // Poll and download all MSAs
    POLL_AND_DOWNLOAD_MSAS(SUBMIT_COLABFOLD.out.tickets)
    
    // Convert to Boltz-2 CSV format
    CONVERT_MSA_TO_BOLTZ_CSV(POLL_AND_DOWNLOAD_MSAS.out.msas)
}



process PREPARE_SEQUENCES {
    output:
    tuple val(params.pred_name), 
          path("sequences.json"), 
          path("protein_chains.txt"), 
          emit: seq_info
    
    script:
    def sequences = []
    
    // Collect sequences that need MSA generation
    ['seq1', 'seq2', 'seq3', 'seq4'].each { prefix ->
        def type = params."${prefix}_type"
        def chain_id = params."${prefix}_chain_id"
        def seq = params."${prefix}_seq"
        def genmsa = params."${prefix}_genmsa"
        
        if (type == 'protein' && seq != null && genmsa == true) {
            sequences.add([
                chain_id: chain_id,
                sequence: seq,
                type: type
            ])
        }
    }
    
    def seq_json = groovy.json.JsonOutput.toJson(sequences)
    def protein_chains = sequences.collect { it.chain_id }.join(',')
    
    """
    echo '${seq_json}' > sequences.json
    echo '${protein_chains}' > protein_chains.txt
    """
}
process SUBMIT_COLABFOLD {
    input:
    tuple val(pred_name), path(seq_json), path(chain_file)
    
    output:
    tuple val(pred_name), path("tickets.json"), path(seq_json), emit: tickets
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import requests
    import sys
    
    # Load sequences
    with open('${seq_json}', 'r') as f:
        sequences = json.load(f)
    
    tickets = []
    
    # Submit each chain separately for unpaired MSAs
    for idx, seq_info in enumerate(sequences):
        chain_id = seq_info['chain_id']
        sequence = seq_info['sequence']
        
        # Boltz uses sequential IDs starting from 101
        fasta_query = f">101\\n{sequence}"
        
        print(f"Submitting unpaired MSA for chain {idx} ({chain_id})", file=sys.stderr)
        
        payload = {
            "q": fasta_query,
            "mode": "env",  # Unpaired mode
        }
        
        response = requests.post(
            'https://api.colabfold.com/ticket/msa',
            data=payload,
            timeout=60
        )
        
        if response.status_code != 200:
            raise Exception(f"API request failed: {response.status_code} - {response.text}")
        
        result = response.json()
        ticket_id = result.get('id')
        
        if not ticket_id:
            raise Exception(f"No ticket ID in response: {result}")
        
        tickets.append({
            'chain_idx': idx,
            'chain_id': chain_id,
            'ticket_id': ticket_id,
            'type': 'unpaired'
        })
        
        print(f"  Ticket ID: {ticket_id}", file=sys.stderr)
    
    # If multiple chains, also submit for paired MSAs
    if len(sequences) > 1:
        print(f"\\nSubmitting paired MSA for all chains", file=sys.stderr)
        
        # Include ALL query sequences in the paired request
        # Boltz assigns sequential IDs starting from 101
        fasta_lines = []
        for idx, seq_info in enumerate(sequences, start=101):
            fasta_lines.append(f">{idx}")
            fasta_lines.append(seq_info['sequence'])
        fasta_query = "\\n".join(fasta_lines)
        
        print(f"Paired query:\\n{fasta_query[:200]}...", file=sys.stderr)
        
        payload = {
            "q": fasta_query,
            "mode": "pairgreedy",  # Paired mode
        }
        
        response = requests.post(
            'https://api.colabfold.com/ticket/pair',  # Note: different endpoint
            data=payload,
            timeout=60
        )
        
        if response.status_code != 200:
            raise Exception(f"Paired API request failed: {response.status_code} - {response.text}")
        
        result = response.json()
        ticket_id = result.get('id')
        
        if not ticket_id:
            raise Exception(f"No ticket ID in paired response: {result}")
        
        tickets.append({
            'chain_idx': -1,  # -1 indicates paired
            'chain_id': 'paired',
            'ticket_id': ticket_id,
            'type': 'paired'
        })
        
        print(f"  Paired Ticket ID: {ticket_id}", file=sys.stderr)
    
    # Save tickets
    with open('tickets.json', 'w') as f:
        json.dump(tickets, f, indent=2)
    
    print(f"\\nTotal tickets: {len(tickets)}", file=sys.stderr)
    """
}


process POLL_AND_DOWNLOAD_MSAS {
    maxRetries 3
    errorStrategy 'retry'
    
    input:
    tuple val(pred_name), path(tickets_json), path(seq_json)
    
    output:
    tuple val(pred_name), 
          path("msas/"),  // Pass the entire directory
          path(seq_json),
          emit: msas
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import time
    from pathlib import Path
    import sys
    import tarfile
    import io
    import json
    
    # Load tickets
    with open('${tickets_json}', 'r') as f:
        tickets = json.load(f)
    
    max_wait = 7200
    poll_interval = 15
    
    for ticket_info in tickets:
        ticket_id = ticket_info['ticket_id']
        ticket_type = ticket_info['type']
        chain_idx = ticket_info['chain_idx']
        
        print(f"\\nProcessing {ticket_type} MSA (ticket: {ticket_id})", file=sys.stderr)
        
        # Poll for completion
        elapsed = 0
        status_url = f'https://api.colabfold.com/ticket/{ticket_id}'
        
        while elapsed < max_wait:
            try:
                response = requests.get(status_url, timeout=30)
                
                if response.status_code != 200:
                    time.sleep(poll_interval)
                    elapsed += poll_interval
                    continue
                
                status_data = response.json()
                status = status_data.get('status', 'UNKNOWN')
                
                print(f"  Status: {status} (elapsed: {elapsed}s)", file=sys.stderr)
                
                if status == 'COMPLETE':
                    break
                elif status == 'ERROR':
                    error_msg = status_data.get('error', 'Unknown error')
                    raise Exception(f'MSA generation failed: {error_msg}')
                elif status in ['PENDING', 'RUNNING']:
                    time.sleep(poll_interval)
                    elapsed += poll_interval
                else:
                    time.sleep(poll_interval)
                    elapsed += poll_interval
                    
            except requests.exceptions.RequestException as e:
                print(f"  Request error: {e}, retrying...", file=sys.stderr)
                time.sleep(poll_interval)
                elapsed += poll_interval
        
        if elapsed >= max_wait:
            raise Exception(f'Timeout waiting for MSA generation')
        
        # Create output directory
        if ticket_type == 'unpaired':
            output_dir = Path(f'msas/chain_{chain_idx}')
        else:
            output_dir = Path('msas/paired')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download results
        download_url = f'https://api.colabfold.com/result/download/{ticket_id}'
        result_response = requests.get(download_url, timeout=120)
        result_response.raise_for_status()
        
        print(f"  Downloaded {len(result_response.content)} bytes", file=sys.stderr)
        
        # Extract tar.gz
        with tarfile.open(fileobj=io.BytesIO(result_response.content), mode='r:gz') as tar:
            for member in tar.getmembers():
                if member.name.endswith('.a3m'):
                    member.name = Path(member.name).name
                    tar.extract(member, path=str(output_dir))
                    print(f"  Extracted: {member.name}", file=sys.stderr)
    
    print(f"\\nAll MSAs downloaded successfully!", file=sys.stderr)
    """
}

process CONVERT_MSA_TO_BOLTZ_CSV {
    publishDir "${params.use_outdir_prefix ? params.outdir_prefix : ''}${params.outdir}/${pred_name}/msas",
               mode: 'copy'
    
    input:
    tuple val(pred_name), path(msas_dir), path(seq_json)
    
    output:
    tuple val(pred_name), path("*.csv"), emit: csv_msas
    
    script:
    """
    #!/usr/bin/env python3
    import json
    from pathlib import Path
    import sys
    import re
    
    # Load sequence info
    with open('${seq_json}', 'r') as f:
        sequences = json.load(f)
    
    num_chains = len(sequences)
    print(f"Processing MSAs for {num_chains} chain(s)", file=sys.stderr)
    
    def parse_a3m_robust(a3m_file):
        sequences_list = []
        
        with open(a3m_file, 'r') as f:
            content = f.read()
        
        # Split by '>' to get each entry
        entries = content.split('>')
        
        for entry in entries:
            if not entry.strip():
                continue
            
            lines = entry.split('\\n', 1)  # Split header from sequence
            if len(lines) < 2:
                continue
                
            header = lines[0].strip()
            sequence_text = lines[1] if len(lines) > 1 else ''
            
            # Remove all newlines and whitespace from sequence
            sequence = ''.join(sequence_text.split())
            
            if sequence:
                sequences_list.append((header, sequence))
        
        return sequences_list
    
    # First, parse the paired MSA structure
    paired_seqs_per_chain = [[] for _ in range(num_chains)]
    
    if num_chains > 1:
        paired_dir = Path('${msas_dir}') / 'paired'
        pair_a3m = paired_dir / 'pair.a3m'
        
        if pair_a3m.exists():
            print(f"\\nParsing paired MSA file...", file=sys.stderr)
            
            # Parse the file
            all_entries = parse_a3m_robust(str(pair_a3m))
            print(f"  Total entries: {len(all_entries)}", file=sys.stderr)
            
            # Find where each chain's section starts
            # Query headers are just numbers: "101", "102", etc.
            chain_starts = []
            for i, (header, seq) in enumerate(all_entries):
                # Check if header is purely numeric (query sequence)
                header_parts = header.split()
                if header_parts and header_parts[0].isdigit():
                    chain_num = int(header_parts[0]) - 101
                    chain_starts.append((chain_num, i))
                    print(f"  Found chain {chain_num} query at index {i} (header: '{header_parts[0]}')", file=sys.stderr)
            
            # Extract paired sequences for each chain
            for idx, (chain_num, start_pos) in enumerate(chain_starts):
                # Find end position (next chain start or end of file)
                if idx + 1 < len(chain_starts):
                    end_pos = chain_starts[idx + 1][1]
                else:
                    end_pos = len(all_entries)
                
                # Extract sequences (skip the query itself at start_pos)
                chain_paired = []
                for i in range(start_pos + 1, end_pos):
                    header, seq = all_entries[i]
                    chain_paired.append((header, seq))
                
                if chain_num < num_chains:
                    paired_seqs_per_chain[chain_num] = chain_paired
                    print(f"  Chain {chain_num}: extracted {len(chain_paired)} paired sequences (indices {start_pos+1} to {end_pos-1})", file=sys.stderr)
    
    # Now process each chain
    for chain_idx in range(num_chains):
        print(f"\\nProcessing chain {chain_idx}", file=sys.stderr)
        
        chain_dir = Path('${msas_dir}') / f'chain_{chain_idx}'
        
        # Load unpaired MSAs for this chain
        unpaired_sequences = []
        
        # Load both uniref and bfd databases
        for db_file in ['uniref.a3m', 'bfd.mgnify30.metaeuk30.smag30.a3m']:
            a3m_path = chain_dir / db_file
            if a3m_path.exists():
                print(f"  Loading {db_file}...", file=sys.stderr)
                seqs = parse_a3m_robust(str(a3m_path))
                print(f"    Found {len(seqs)} sequences", file=sys.stderr)
                unpaired_sequences.extend(seqs)
            else:
                print(f"  Warning: {db_file} not found", file=sys.stderr)
        
        # Get paired sequences for this chain
        paired_sequences = paired_seqs_per_chain[chain_idx]
        
        # Build CSV for this chain
        # Order: Query (key=0), Paired (keys=1,2,3...), Unpaired (keys=-1)
        output_csv = f"${pred_name}_{chain_idx}.csv"
        
        with open(output_csv, 'w') as csvfile:
            # Write header
            csvfile.write('key,sequence\\n')
            
            # 1. Query sequence (key = 0)
            if unpaired_sequences:
                query_header, query_seq = unpaired_sequences[0]
                csvfile.write(f'0,{query_seq}\\n')
                print(f"  Query ({len(query_seq)} aa): {query_seq[:50]}...", file=sys.stderr)
            
            # 2. Paired sequences (keys = 1, 2, 3, ...)
            for pair_idx, (header, seq) in enumerate(paired_sequences, start=1):
                csvfile.write(f'{pair_idx},{seq}\\n')
            
            # 3. Unpaired sequences (all key = -1)
            # Skip the first one (query)
            for header, seq in unpaired_sequences[1:]:
                csvfile.write(f'-1,{seq}\\n')
        
        paired_count = len(paired_sequences)
        unpaired_count = len(unpaired_sequences) - 1  # Subtract query
        total = 1 + paired_count + unpaired_count
        
        print(f"  Created {output_csv}:", file=sys.stderr)
        print(f"    Query (key=0): 1", file=sys.stderr)
        print(f"    Paired (keys=1-{paired_count}): {paired_count}", file=sys.stderr)
        print(f"    Unpaired (key=-1): {unpaired_count}", file=sys.stderr)
        print(f"    Total: {total}", file=sys.stderr)
    
    # Verify outputs
    csv_files = list(Path('.').glob('*.csv'))
    print(f"\\nCreated {len(csv_files)} CSV files", file=sys.stderr)
    for csv_file in csv_files:
        line_count = sum(1 for _ in open(csv_file)) - 1
        size_kb = csv_file.stat().st_size / 1024
        print(f"  {csv_file.name} ({line_count} sequences, {size_kb:.1f} KB)", file=sys.stderr)
    """
}









