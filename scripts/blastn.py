#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import subprocess
import log,traceback


def _run(cmd):
    """Run a command list, raising on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            'Command failed (exit %d): %s\nstderr: %s'
            % (result.returncode, ' '.join(cmd), result.stderr)
        )


def blastn(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        _run(['blastn', '-db', db_path, '-query', q_path,
              '-evalue', str(params.blastn_evalue),
              '-perc_identity', str(params.blastn_ident),
              '-word_size', str(params.blastn_word_size),
              '-num_threads', str(args.p),
              '-outfmt', '6', '-out', outfpath])
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def blastn_culling(args, params, q_path, db_path, outfpath, culling_limit):
    log.logger.debug('started')
    try:
        _run(['blastn', '-db', db_path, '-query', q_path,
              '-evalue', str(params.blastn_evalue),
              '-perc_identity', str(params.blastn_ident),
              '-word_size', str(params.blastn_word_size),
              '-num_threads', str(args.p),
              '-outfmt', '6', '-out', outfpath,
              '-culling_limit', str(culling_limit)])
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def blastn_single_thread(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        _run(['blastn', '-db', db_path, '-query', q_path,
              '-evalue', str(params.blastn_evalue),
              '-perc_identity', str(params.blastn_ident),
              '-word_size', str(params.blastn_word_size),
              '-num_threads', '1',
              '-outfmt', '6', '-out', outfpath])
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def blastn_for_unknown_rep_ident(args, params, q_path, db_path, outfpath):
    log.logger.debug('started')
    try:
        _run(['blastn', '-db', db_path, '-query', q_path,
              '-evalue', str(params.blastn_evalue),
              '-perc_identity', str(params.blastn_ident),
              '-word_size', '7',
              '-num_threads', str(args.p),
              '-outfmt', '6', '-out', outfpath])
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def makeblastdb(fasta_file, dbpath):
    log.logger.debug('started')
    try:
        _run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl',
              '-out', dbpath, '-parse_seqids'])
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
