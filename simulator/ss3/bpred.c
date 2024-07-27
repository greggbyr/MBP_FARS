/* bpred.c - branch predictor routines */

/* SimpleScalar(TM) Tool Suite
 * Copyright (C) 1994-2003 by Todd M. Austin, Ph.D. and SimpleScalar, LLC.
 * All Rights Reserved. 
 * 
 * THIS IS A LEGAL DOCUMENT, BY USING SIMPLESCALAR,
 * YOU ARE AGREEING TO THESE TERMS AND CONDITIONS.
 * 
 * No portion of this work may be used by any commercial entity, or for any
 * commercial purpose, without the prior, written permission of SimpleScalar,
 * LLC (info@simplescalar.com). Nonprofit and noncommercial use is permitted
 * as described below.
 * 
 * 1. SimpleScalar is provided AS IS, with no warranty of any kind, express
 * or implied. The user of the program accepts full responsibility for the
 * application of the program and the use of any results.
 * 
 * 2. Nonprofit and noncommercial use is encouraged. SimpleScalar may be
 * downloaded, compiled, executed, copied, and modified solely for nonprofit,
 * educational, noncommercial research, and noncommercial scholarship
 * purposes provided that this notice in its entirety accompanies all copies.
 * Copies of the modified software can be delivered to persons who use it
 * solely for nonprofit, educational, noncommercial research, and
 * noncommercial scholarship purposes provided that this notice in its
 * entirety accompanies all copies.
 * 
 * 3. ALL COMMERCIAL USE, AND ALL USE BY FOR PROFIT ENTITIES, IS EXPRESSLY
 * PROHIBITED WITHOUT A LICENSE FROM SIMPLESCALAR, LLC (info@simplescalar.com).
 * 
 * 4. No nonprofit user may place any restrictions on the use of this software,
 * including as modified by the user, by any other authorized user.
 * 
 * 5. Noncommercial and nonprofit users may distribute copies of SimpleScalar
 * in compiled or executable form as set forth in Section 2, provided that
 * either: (A) it is accompanied by the corresponding machine-readable source
 * code, or (B) it is accompanied by a written offer, with no time limit, to
 * give anyone a machine-readable copy of the corresponding source code in
 * return for reimbursement of the cost of distribution. This written offer
 * must permit verbatim duplication by anyone, or (C) it is distributed by
 * someone who received only the executable form, and is accompanied by a
 * copy of the written offer of source code.
 * 
 * 6. SimpleScalar was developed by Todd M. Austin, Ph.D. The tool suite is
 * currently maintained by SimpleScalar LLC (info@simplescalar.com). US Mail:
 * 2395 Timbercrest Court, Ann Arbor, MI 48105.
 * 
 * Copyright (C) 1994-2003 by Todd M. Austin, Ph.D. and SimpleScalar, LLC.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "host.h"
#include "misc.h"
#include "machine.h"
#include "bpred.h"

/* turn this on to enable the SimpleScalar 2.0 RAS bug */
/* #define RAS_BUG_COMPATIBLE */

/* create a branch predictor */
struct bpred_t *			/* branch predictory instance */
bpred_create(enum bpred_class class,	/* type of predictor to create */
	     unsigned int bimod_size,	/* bimod table size */
	     unsigned int l1size,	/* 2lev l1 table size */
	     unsigned int l2size,	/* 2lev l2 table size */
	     unsigned int meta_size,	/* meta table size */
	     unsigned int shift_width,	/* history register width */
	     unsigned int xor,  	/* history xor address flag */
     	     unsigned int head_table_width, /*TS head table width*/
	     unsigned int cht_size, /*CHBP correctness history table width*/
	     unsigned int btb_sets,	/* number of sets in BTB */ 
	     unsigned int btb_assoc,	/* BTB associativity */
	     unsigned int retstack_size) /* num entries in ret-addr stack */
	     //unsigned int ts_enabled)	/* TSBP Enabled Flag */
{
  struct bpred_t *pred;

  if (!(pred = calloc(1, sizeof(struct bpred_t))))
    fatal("out of virtual memory");

  pred->class = class;

  switch (class) {
  case BPredComb:
    /* bimodal component */
    pred->fwd_dirpred.bimod = 
      bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);

    //Need reverse bimod now too
    pred->rev_dirpred.bimod =
      bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);

    /* 2-level component */
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	pred->rev_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);

	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(pred, (shift_width + 1));

    /* metapredictor component */
    pred->fwd_dirpred.meta = 
      bpred_dir_create(BPred2bit, meta_size, 0, 0, 0);
	  
	//Need reverse meta too
	pred->rev_dirpred.meta = 
      bpred_dir_create(BPred2bit, meta_size, 0, 0, 0);

    break;

  case BPred2Level:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(class, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	pred->rev_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(pred, (shift_width + 1));
	
    break;
	
  case BPredOB:
	pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	pred->rev_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(pred, (shift_width + 1));
	
	// Outcome Buffer (OB)
	bpred_ob_create(pred, 16384); /* Hard setting OB width to 16k for now*/
	
    break;
	
  case BPredOHT:
	pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	pred->rev_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(pred, (shift_width + 1));
	
	//Need FWD and REV OHTs
	pred->fwd_dirpred.oht = 
      bpred_oht_create(class, l2size);
	  
	pred->rev_dirpred.oht = 
      bpred_oht_create(class, l2size);
	
    break;
	
  case BPredMBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	pred->rev_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(pred, (shift_width + 1));
	
	//Need FWD and REV OHTs
	pred->fwd_dirpred.oht = 
      bpred_oht_create(class, l2size);
	  
	pred->rev_dirpred.oht = 
      bpred_oht_create(class, l2size);
	
	// Outcome Buffer (OB)
	bpred_ob_create(pred, 16384); /* Hard setting OB width to 16k for now*/
	
    break;

  case BPredTSBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	
	pred->fwd_dirpred.tsbp =
	  bpred_ts_create(class, 1, head_table_width, ((unsigned int)l2size << 3)); /* md_addr_t is the size of the PC*/
    break;
	
  case BPredCHBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	
	pred->fwd_dirpred.chbp =
	  bpred_chbp_create(class, 1, cht_size);
    break;

  case BPred2bit:
    pred->fwd_dirpred.bimod = 
      bpred_dir_create(class, bimod_size, 0, 0, 0);

    //Need reverse bimod now too
    pred->rev_dirpred.bimod =
      bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);
  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus predictor class");
  }

  /* allocate ret-addr stack */
  switch (class) {
  case BPredComb:
  case BPred2Level:
  case BPredTSBP:
  case BPredCHBP:
  case BPred2bit:
  case BPredOB:
  case BPredOHT:
  case BPredMBP:
    {
      int i;

      /* allocate BTB */
      if (!btb_sets || (btb_sets & (btb_sets-1)) != 0)
	fatal("number of BTB sets must be non-zero and a power of two");
      if (!btb_assoc || (btb_assoc & (btb_assoc-1)) != 0)
	fatal("BTB associativity must be non-zero and a power of two");

      if (!(pred->btb.btb_data = calloc(btb_sets * btb_assoc,
					sizeof(struct bpred_btb_ent_t))))
	fatal("cannot allocate BTB");

      pred->btb.sets = btb_sets;
      pred->btb.assoc = btb_assoc;

      if (pred->btb.assoc > 1)
	for (i=0; i < (pred->btb.assoc*pred->btb.sets); i++)
	  {
	    if (i % pred->btb.assoc != pred->btb.assoc - 1)
	      pred->btb.btb_data[i].next = &pred->btb.btb_data[i+1];
	    else
	      pred->btb.btb_data[i].next = NULL;
	    
	    if (i % pred->btb.assoc != pred->btb.assoc - 1)
	      pred->btb.btb_data[i+1].prev = &pred->btb.btb_data[i];
	  }

      /* allocate retstack */
      if ((retstack_size & (retstack_size-1)) != 0)
	fatal("Return-address-stack size must be zero or a power of two");
      
      pred->retstack.size = retstack_size;
      if (retstack_size)
	if (!(pred->retstack.stack = calloc(retstack_size, 
					    sizeof(struct bpred_btb_ent_t))))
	  fatal("cannot allocate return-address-stack");
      pred->retstack.tos = retstack_size - 1;
      
      break;
    }

  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus predictor class");
  }

  return pred;
}

/* Create an instance of the Outcome Buffer (OB) */
void bpred_ob_create (struct bpred_t *pred, int width) {
	// Outcome Buffer (OB)
	pred->ob.width = width; /* hardcoding to 16 kb for now, will make configurable later */
	
	pred->ob.end = width - 1;	/* End of buffer window at width - 1 */
	
	pred->ob.beg = 0;	/* Beginning of buffer window at 0 */
	
	if (!(pred->ob.fv = calloc(width, sizeof(int))))
		fatal("cannot allocate OB FV bits");
	
	if (!(pred->ob.rv = calloc(width, sizeof(int))))
		fatal("cannot allocate OB RV bits");
	
	if (!(pred->ob.oc = calloc(width, sizeof(int))))
		fatal("cannot allocate OB OC bits");
}

/* Create an instance of the Future History Buffer (FHB) */
void bpred_fhb_create (struct bpred_t *pred, int size) {
	pred->fhb.size = size;
	
	pred->fhb.top = size - 1;
	
	pred->fhb.bot = 0;
	
	if (!(pred->fhb.fv = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB FV bits");
	
	if (!(pred->fhb.rv = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB RV bits");
	
	if (!(pred->fhb.o = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB O bits");
	
	if (!(pred->fhb.addr = calloc(size, sizeof(md_addr_t))))
		fatal("cannot allocate FHB ADDR bits");
}

/* create a branch direction predictor */
struct bpred_dir_t *		/* branch direction predictor instance */
bpred_dir_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int l1size,	 	/* level-1 table size */
  unsigned int l2size,	 	/* level-2 table size (if relevant) */
  unsigned int shift_width,	/* history register width */
  unsigned int xor)	    	/* history xor address flag */
{
  struct bpred_dir_t *pred_dir;
  unsigned int cnt;
  int flipflop;

  if (!(pred_dir = calloc(1, sizeof(struct bpred_dir_t))))
    fatal("out of virtual memory");

  pred_dir->class = class;

  cnt = -1;
  switch (class) {
  case BPred2Level:
    {
      if (!l1size || (l1size & (l1size-1)) != 0)
	fatal("level-1 size, `%d', must be non-zero and a power of two", 
	      l1size);
      pred_dir->config.two.l1size = l1size;
      
      if (!l2size || (l2size & (l2size-1)) != 0)
	fatal("level-2 size, `%d', must be non-zero and a power of two", 
	      l2size);
      pred_dir->config.two.l2size = l2size;
      
      if (!shift_width || shift_width > 30)
	fatal("shift register width, `%d', must be non-zero and positive",
	      shift_width);
      pred_dir->config.two.shift_width = shift_width;
      
      pred_dir->config.two.xor = xor;
      pred_dir->config.two.shiftregs = calloc(l1size, sizeof(int));
      if (!pred_dir->config.two.shiftregs)
	fatal("cannot allocate shift register table");
      
      pred_dir->config.two.l2table = calloc(l2size, sizeof(unsigned char));
      if (!pred_dir->config.two.l2table)
	fatal("cannot allocate second level table");

      /* initialize counters to weakly this-or-that */
      flipflop = 1;
      for (cnt = 0; cnt < l2size; cnt++)
	{
	  pred_dir->config.two.l2table[cnt] = flipflop;
	  flipflop = 3 - flipflop;
	}

      break;
    }

  case BPred2bit:
    if (!l1size || (l1size & (l1size-1)) != 0)
      fatal("2bit table size, `%d', must be non-zero and a power of two", 
	    l1size);
    pred_dir->config.bimod.size = l1size;
    if (!(pred_dir->config.bimod.table =
	  calloc(l1size, sizeof(unsigned char))))
      fatal("cannot allocate 2bit storage");
    /* initialize counters to weakly this-or-that */
    flipflop = 1;
    for (cnt = 0; cnt < l1size; cnt++)
      {
	pred_dir->config.bimod.table[cnt] = flipflop;
	flipflop = 3 - flipflop;
      }

    break;

  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus branch direction predictor class");
  }

  return pred_dir;
}

/* create a branch direction predictor */
struct bpred_ts_t *		/* temporal stream branch predictor instance */
bpred_ts_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int ts_enabled,              /* TSBP Enabled Flag */
  unsigned int head_table_width,	 	/* head table width */
  unsigned int head_table_size)			/* head table size */
{
  struct bpred_ts_t *pred_ts;
  unsigned int key;
  
  if (!(pred_ts = calloc(1, sizeof(struct bpred_ts_t))))
    fatal("out of virtual memory");

  pred_ts->class = class;

  if (!head_table_size || (head_table_size & (head_table_size - 1)) != 0)
	fatal("head table size, `%d', must be non-zero and a power of two", head_table_size);
  
  pred_ts->ts.head_table_size = head_table_size;
  
  if (!head_table_width || head_table_width > 30)
	fatal("head table width, `%d', must be non-zero and positive", head_table_width);
  
  pred_ts->ts.head_table_width = head_table_width;
 
  pred_ts->ts.head_table = calloc(head_table_size, sizeof(unsigned int));
  //if (head_table_width % 8) { 
  //	pred_ts->ts.head_table = calloc(head_table_size, ((head_table_width / 8) + 1));
  //} else {
//	pred_ts->ts.head_table = calloc(head_table_size, (head_table_width / 8));
  //}
  
  if (!pred_ts->ts.head_table)
	fatal("cannot allocate head table");

  /* initializing head table entries to NULL*/
  for (key = 0; key < head_table_size; key++)
	  pred_ts->ts.head_table[key] = NULL;

  pred_ts->ts.correctness_width = 2 << (head_table_width - 1);
  
  pred_ts->ts.correctness_buffer = calloc(pred_ts->ts.correctness_width, sizeof(bool_t));

  if (!pred_ts->ts.correctness_buffer)
	fatal("cannot allocate correctness buffer");

   /* initializing CB bits to 1*/
  for (key = 0; key < pred_ts->ts.correctness_width; key++)
  	pred_ts->ts.correctness_buffer[key] = 1;

  /* initialize current head of the correctness buffer */
  pred_ts->ts.head = 0;

  /* initialize current tail of the correctness buffer */
  pred_ts->ts.tail = 0;
  
  /* initialize replay to false */
  pred_ts->ts.replay = FALSE;

  /* initialize enabled flag */
  if (ts_enabled == 0)
          pred_ts->ts.enabled = FALSE;
  else
          pred_ts->ts.enabled = TRUE;

  return pred_ts;
}

/* create a branch direction predictor */
struct bpred_chbp_t *		/* temporal stream branch predictor instance */
bpred_chbp_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int chbp_enabled,              /* CHBP Enabled Flag */
  unsigned int cht_size)			/* Correctness History Table size */
{
  struct bpred_chbp_t *pred_chbp;
  //unsigned int key;
  
  if (!(pred_chbp = calloc(1, sizeof(struct bpred_chbp_t))))
    fatal("out of virtual memory");

  pred_chbp->class = class;

  if (!cht_size || (cht_size & (cht_size - 1)) != 0)
	fatal("correctness history table size, `%d', must be non-zero and a power of two", cht_size);
  
  pred_chbp->chbp.cht_size = cht_size;
  
  /* Need to allocate for source PC, replay bit, correctness bit, valid bit, and destination PC*/
  pred_chbp->chbp.cht_spc = calloc(cht_size, sizeof(md_addr_t)); 
  
  if (!pred_chbp->chbp.cht_spc)
	fatal("cannot allocate correctness history table source pc bits");

  pred_chbp->chbp.cht_replay = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_replay)
	fatal("cannot allocate correctness history table replay bits");

  pred_chbp->chbp.cht_correct = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_correct)
	fatal("cannot allocate correctness history table correct bits");

  pred_chbp->chbp.cht_valid = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_valid)
	fatal("cannot allocate correctness history table valid bits");

  pred_chbp->chbp.cht_dpc = calloc(cht_size, sizeof(md_addr_t)); 
  
  if (!pred_chbp->chbp.cht_dpc)
	fatal("cannot allocate correctness history table destination pc bits");

  /* initialize enabled flag */
  if (chbp_enabled == 0)
          pred_chbp->chbp.enabled = FALSE;
  else
          pred_chbp->chbp.enabled = TRUE;

  return pred_chbp;
}

/* create an outcome history table */
struct bpred_oht_t *		/* outcome history table instance */
bpred_oht_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int oht_size)			/* Outcome History Table size */
{
  struct bpred_oht_t *pred_oht;
  //unsigned int key;
  
  if (!(pred_oht = calloc(1, sizeof(struct bpred_oht_t))))
    fatal("out of virtual memory");

  pred_oht->class = class;

  if (!oht_size || (oht_size & (oht_size - 1)) != 0)
	fatal("outcome history table size, `%d', must be non-zero and a power of two", oht_size);
  
  pred_oht->oht.size = oht_size;
  
  /* Need to allocate for outcome (OC) and valid bits */
  pred_oht->oht.oc = calloc(oht_size, sizeof(int)); 
  
  if (!pred_oht->oht.oc)
	fatal("cannot allocate outcome history table outcome (pc) bits");

  pred_oht->oht.valid = calloc(oht_size, sizeof(int)); 
  
  if (!pred_oht->oht.valid)
	fatal("cannot allocate outcome history table valid bits");

  return pred_oht;
}

/* print branch direction predictor configuration */
void
bpred_dir_config(
  struct bpred_dir_t *pred_dir,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  switch (pred_dir->class) {
  case BPred2Level:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredOB:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredOHT:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredMBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredTSBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;
	
  case BPredCHBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPred2bit:
    fprintf(stream, "pred_dir: %s: 2-bit: %d entries, direct-mapped\n",
      name, pred_dir->config.bimod.size);
    break;

  case BPredTaken:
    fprintf(stream, "pred_dir: %s: predict taken\n", name);
    break;

  case BPredNotTaken:
    fprintf(stream, "pred_dir: %s: predict not taken\n", name);
    break;

  default:
    panic("bogus branch direction predictor class");
  }
}

/* print temporal stream predictor configuration */
void
bpred_ts_config(
  struct bpred_ts_t *pred_ts,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  fprintf(stream,
    "pred_ts: %s: tsbp: %d ht-sz, %d ht-wd, %d cr-wd\n",
    name, pred_ts->ts.head_table_size, pred_ts->ts.head_table_width,
	pred_ts->ts.correctness_width);
}

/* print mississippi branch predictor configuration */
void
bpred_chbp_config(
  struct bpred_chbp_t *pred_chbp,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  fprintf(stream,
    "pred_chbp: %s: chbp: %d cht-sz\n",
    name, pred_chbp->chbp.cht_size);
}

/* print branch predictor configuration */
void
bpred_config(struct bpred_t *pred,	/* branch predictor instance */
	     FILE *stream)		/* output stream */
{
  switch (pred->class) {
  case BPredComb:
    bpred_dir_config (pred->fwd_dirpred.bimod, "bimod", stream);
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    bpred_dir_config (pred->fwd_dirpred.meta, "meta", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPred2Level:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
	
  case BPredOB:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredOHT:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
  
  case BPredMBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredTSBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
	bpred_ts_config (pred->fwd_dirpred.tsbp, "tsbp", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
	
  case BPredCHBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
	bpred_chbp_config (pred->fwd_dirpred.chbp, "chbp", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPred2bit:
    bpred_dir_config (pred->fwd_dirpred.bimod, "bimod", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredTaken:
    bpred_dir_config (pred->fwd_dirpred.bimod, "taken", stream);
    break;
  case BPredNotTaken:
    bpred_dir_config (pred->fwd_dirpred.bimod, "nottaken", stream);
    break;

  default:
    panic("bogus branch predictor class");
  }
}

/* print predictor stats */
void
bpred_stats(struct bpred_t *pred,	/* branch predictor instance */
	    FILE *stream)		/* output stream */
{
	/* stats */
  fprintf(stream, "pred: addr-prediction rate = %f\n",
	  (double)pred->addr_hits/(double)(pred->addr_hits+pred->misses));
  fprintf(stream, "pred: dir-prediction rate = %f\n",
	  (double)pred->dir_hits/(double)(pred->dir_hits+pred->misses));
  
	/* reverse stats */
  fprintf(stream, "pred: reverse addr-prediction rate = %f\n",
	  (double)pred->reverse_addr_hits/(double)(pred->reverse_addr_hits+pred->reverse_misses));
  fprintf(stream, "pred: reverse dir-prediction rate = %f\n",
	  (double)pred->reverse_dir_hits/(double)(pred->reverse_dir_hits+pred->reverse_misses));
}

/* register branch predictor stats */
void
bpred_reg_stats(struct bpred_t *pred,	/* branch predictor instance */
		struct stat_sdb_t *sdb)	/* stats database */
{
  char buf[512], buf1[512], *name;

  /* get a name for this predictor */
  switch (pred->class)
    {
    case BPredComb:
      name = "bpred_comb";
      break;
    case BPred2Level:
      name = "bpred_2lev";
      break;
     case BPredOB:
      name = "bpred_ob";
      break;
     case BPredOHT:
      name = "bpred_oht";
      break;
     case BPredMBP:
      name = "bpred_mbp";
      break;
	case BPredTSBP:
      name = "bpred_tsbp";
      break;
	case BPredCHBP:
      name = "bpred_chbp";
      break;
    case BPred2bit:
      name = "bpred_bimod";
      break;
    case BPredTaken:
      name = "bpred_taken";
      break;
    case BPredNotTaken:
      name = "bpred_nottaken";
      break;
    default:
      panic("bogus branch predictor class");
    }

  sprintf(buf, "%s.lookups", name);
  stat_reg_counter(sdb, buf, "total number of bpred lookups",
		   &pred->lookups, 0, NULL);
  sprintf(buf, "%s.updates", name);
  sprintf(buf1, "%s.dir_hits + %s.misses", name, name);
  stat_reg_formula(sdb, buf, "total number of updates", buf1, "%12.0f");
  sprintf(buf, "%s.addr_hits", name);
  stat_reg_counter(sdb, buf, "total number of address-predicted hits", 
		   &pred->addr_hits, 0, NULL);
  sprintf(buf, "%s.dir_hits", name);
  stat_reg_counter(sdb, buf, 
		   "total number of direction-predicted hits "
		   "(includes addr-hits)", 
		   &pred->dir_hits, 0, NULL);
  if (pred->class == BPredComb)
    {
      sprintf(buf, "%s.used_bimod", name);
      stat_reg_counter(sdb, buf, 
		       "total number of bimodal predictions used", 
		       &pred->used_bimod, 0, NULL);
      sprintf(buf, "%s.used_2lev", name);
      stat_reg_counter(sdb, buf, 
		       "total number of 2-level predictions used", 
		       &pred->used_2lev, 0, NULL);
    }
  if (pred->class == BPredTSBP || pred->class == BPredCHBP) {
       sprintf(buf, "%s.replays", name);
       stat_reg_counter(sdb, buf, "total number of replays", &pred->replays, 0, NULL);
  }
  sprintf(buf, "%s.misses", name);
  stat_reg_counter(sdb, buf, "total number of misses", &pred->misses, 0, NULL);
  sprintf(buf, "%s.jr_hits", name);
  stat_reg_counter(sdb, buf,
		   "total number of address-predicted hits for JR's",
		   &pred->jr_hits, 0, NULL);
  sprintf(buf, "%s.jr_seen", name);
  stat_reg_counter(sdb, buf,
		   "total number of JR's seen",
		   &pred->jr_seen, 0, NULL);
  sprintf(buf, "%s.jr_non_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of address-predicted hits for non-RAS JR's",
		   &pred->jr_non_ras_hits, 0, NULL);
  sprintf(buf, "%s.jr_non_ras_seen.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of non-RAS JR's seen",
		   &pred->jr_non_ras_seen, 0, NULL);
  sprintf(buf, "%s.bpred_addr_rate", name);
  sprintf(buf1, "%s.addr_hits / %s.updates", name, name);
  stat_reg_formula(sdb, buf,
		   "branch address-prediction rate (i.e., addr-hits/updates)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.bpred_dir_rate", name);
  sprintf(buf1, "%s.dir_hits / %s.updates", name, name);
  stat_reg_formula(sdb, buf,
		  "branch direction-prediction rate (i.e., all-hits/updates)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.bpred_jr_rate", name);
  sprintf(buf1, "%s.jr_hits / %s.jr_seen", name, name);
  stat_reg_formula(sdb, buf,
		  "JR address-prediction rate (i.e., JR addr-hits/JRs seen)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.bpred_jr_non_ras_rate.PP", name);
  sprintf(buf1, "%s.jr_non_ras_hits.PP / %s.jr_non_ras_seen.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "non-RAS JR addr-pred rate (ie, non-RAS JR hits/JRs seen)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.retstack_pushes", name);
  stat_reg_counter(sdb, buf,
		   "total number of address pushed onto ret-addr stack",
		   &pred->retstack_pushes, 0, NULL);
  sprintf(buf, "%s.retstack_pops", name);
  stat_reg_counter(sdb, buf,
		   "total number of address popped off of ret-addr stack",
		   &pred->retstack_pops, 0, NULL);
  sprintf(buf, "%s.used_ras.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of RAS predictions used",
		   &pred->used_ras, 0, NULL);
  sprintf(buf, "%s.ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of RAS hits",
		   &pred->ras_hits, 0, NULL);
  sprintf(buf, "%s.ras_rate.PP", name);
  sprintf(buf1, "%s.ras_hits.PP / %s.used_ras.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "RAS prediction rate (i.e., RAS hits/used RAS)",
		   buf1, "%9.4f");
		   
		   
		  /*  Reverse Data Sets */
		   
  sprintf(buf, "%s.reverse_lookups", name);
  stat_reg_counter(sdb, buf, "total number of reverse bpred lookups",
		   &pred->reverse_lookups, 0, NULL);
  sprintf(buf, "%s.reverse_updates", name);
  sprintf(buf1, "%s.reverse_dir_hits + %s.reverse_misses", name, name);
  stat_reg_formula(sdb, buf, "total number of reverse updates", buf1, "%12.0f");
  sprintf(buf, "%s.reverse_addr_hits", name);
  stat_reg_counter(sdb, buf, "total number of reverse address-predicted hits", 
		   &pred->reverse_addr_hits, 0, NULL);
  sprintf(buf, "%s.reverse_dir_hits", name);
  stat_reg_counter(sdb, buf, 
		   "total number of reverse direction-predicted hits "
		   "(includes addr-hits)", 
		   &pred->reverse_dir_hits, 0, NULL);
  if (pred->class == BPredComb)
    {
      sprintf(buf, "%s.reverse_used_bimod", name);
      stat_reg_counter(sdb, buf, 
		       "total number of reverse bimodal predictions used", 
		       &pred->reverse_used_bimod, 0, NULL);
      sprintf(buf, "%s.reverse_used_2lev", name);
      stat_reg_counter(sdb, buf, 
		       "total number of reverse 2-level predictions used", 
		       &pred->reverse_used_2lev, 0, NULL);
    }
  if (pred->class == BPredTSBP || pred->class == BPredCHBP) {
       sprintf(buf, "%s.reverse_replays", name);
       stat_reg_counter(sdb, buf, "total number of reverse replays", &pred->reverse_replays, 0, NULL);
  }
  sprintf(buf, "%s.reverse_misses", name);
  stat_reg_counter(sdb, buf, "total number of reverse misses", &pred->reverse_misses, 0, NULL);
  sprintf(buf, "%s.reverse_jr_hits", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address-predicted hits for JR's",
		   &pred->reverse_jr_hits, 0, NULL);
  sprintf(buf, "%s.reverse_jr_seen", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse JR's seen",
		   &pred->reverse_jr_seen, 0, NULL);
  sprintf(buf, "%s.reverse_jr_non_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address-predicted hits for non-RAS JR's",
		   &pred->reverse_jr_non_ras_hits, 0, NULL);
  sprintf(buf, "%s.reverse_jr_non_ras_seen.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse non-RAS JR's seen",
		   &pred->reverse_jr_non_ras_seen, 0, NULL);
  sprintf(buf, "%s.reverse_bpred_addr_rate", name);
  sprintf(buf1, "%s.reverse_addr_hits / %s.reverse_updates", name, name);
  stat_reg_formula(sdb, buf,
		   "branch address-prediction rate (i.e., addr-hits/updates)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_dir_rate", name);
  sprintf(buf1, "%s.reverse_dir_hits / %s.reverse_updates", name, name);
  stat_reg_formula(sdb, buf,
		  "branch direction-prediction rate (i.e., all-hits/updates)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_jr_rate", name);
  sprintf(buf1, "%s.reverse_jr_hits / %s.reverse_jr_seen", name, name);
  stat_reg_formula(sdb, buf,
		  "JR address-prediction rate (i.e., JR addr-hits/JRs seen)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_jr_non_ras_rate.PP", name);
  sprintf(buf1, "%s.reverse_jr_non_ras_hits.PP / %s.reverse_jr_non_ras_seen.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "non-RAS JR addr-pred rate (ie, non-RAS JR hits/JRs seen)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.reverse_retstack_pushes", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address pushed onto ret-addr stack",
		   &pred->reverse_retstack_pushes, 0, NULL);
  sprintf(buf, "%s.reverse_retstack_pops", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address popped off of ret-addr stack",
		   &pred->reverse_retstack_pops, 0, NULL);
  sprintf(buf, "%s.reverse_used_ras.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse RAS predictions used",
		   &pred->reverse_used_ras, 0, NULL);
  sprintf(buf, "%s.reverse_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse RAS hits",
		   &pred->reverse_ras_hits, 0, NULL);
  sprintf(buf, "%s.reverse_ras_rate.PP", name);
  sprintf(buf1, "%s.reverse_ras_hits.PP / %s.reverse_used_ras.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "RAS prediction rate (i.e., RAS hits/used RAS)",
		   buf1, "%9.4f");
}

void
bpred_after_priming(struct bpred_t *bpred)
{
  if (bpred == NULL)
    return;

	/* Stats */
  bpred->lookups = 0;
  bpred->addr_hits = 0;
  bpred->dir_hits = 0;
  bpred->used_ras = 0;
  bpred->used_bimod = 0;
  bpred->used_2lev = 0;
  bpred->jr_hits = 0;
  bpred->jr_seen = 0;
  bpred->misses = 0;
  bpred->replays = 0;
  bpred->retstack_pops = 0;
  bpred->retstack_pushes = 0;
  bpred->ras_hits = 0;
  
	/* Reverse Stats */ //Set to 1 for now to avoid runtime issues
  bpred->reverse_lookups = 0;
  bpred->reverse_addr_hits = 0;
  bpred->reverse_dir_hits = 0;
  bpred->reverse_used_ras = 0;
  bpred->reverse_used_bimod = 0;
  bpred->reverse_used_2lev = 0;
  bpred->reverse_jr_hits = 0;
  bpred->reverse_jr_seen = 0;
  bpred->reverse_misses = 0;
  bpred->reverse_replays = 0;
  bpred->reverse_retstack_pops = 0;
  bpred->reverse_retstack_pushes = 0;
  bpred->reverse_ras_hits = 0;
}

#define BIMOD_HASH(PRED, ADDR)						\
  ((((ADDR) >> 19) ^ ((ADDR) >> MD_BR_SHIFT)) & ((PRED)->config.bimod.size-1))
    /* was: ((baddr >> 16) ^ baddr) & (pred->fwd_dirpred.bimod.size-1) */

/* Used to calculate 2nd level table indexes/keys for 2lev, tsbp, and chbp*/

int key_from_features (struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
		 md_addr_t baddr)		/* branch address */
{
	int l1index, key;

    /* traverse 2-level tables */
    l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
    key = pred_dir->config.two.shiftregs[l1index];
        
	if (pred_dir->config.two.xor) {
	    /* this L2 index computation is more "compatible" to McFarling's
	       verison of it, i.e., if the PC xor address component is only
	       part of the index, take the lower order address bits for the
	       other part of the index, rather than the higher order ones */
#if 1  
	    key = (((key ^ (baddr >> MD_BR_SHIFT))
			& ((1 << pred_dir->config.two.shift_width) - 1))
		    | ((baddr >> MD_BR_SHIFT)
			<< pred_dir->config.two.shift_width));
#else   
		key = key ^ (baddr >> MD_BR_SHIFT);
#endif
	} else {
		key = key | ((baddr >> MD_BR_SHIFT) << pred_dir->config.two.shift_width);
	}
	
	return key;
}

/* predicts a branch direction */
char *						/* pointer to counter */
bpred_dir_lookup(struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
		 md_addr_t baddr)		/* branch address */
{
  unsigned char *p = NULL;

  /* Except for jumps, get a pointer to direction-prediction bits */
  switch (pred_dir->class) {
    case BPred2Level:
    case BPredOB:
    case BPredOHT:
    case BPredMBP:  
    case BPredTSBP:         /*Add TSBP case, should be same as 2 level to get base prediction*/
	case BPredCHBP:         /*Add CHBP case, should be same as 2 level to get base prediction*/
      {
		int l2index = key_from_features (pred_dir, baddr); //Get un-masked l2index  
		
        l2index = l2index & (pred_dir->config.two.l2size - 1);

        /* get a pointer to prediction state information */
        p = &pred_dir->config.two.l2table[l2index];
      }
      break;
    case BPred2bit:
      p = &pred_dir->config.bimod.table[BIMOD_HASH(pred_dir, baddr)];
      break;
    case BPredTaken:
    case BPredNotTaken:
      break;
    default:
      panic("bogus branch direction predictor class");
    }

  return (char *)p;
}

/* probe a predictor for a next fetch address, the predictor is probed
   with branch address BADDR, the branch target is BTARGET (used for
   static predictors), and OP is the instruction opcode (used to simulate
   predecode bits; a pointer to the predictor state entry (or null for jumps)
   is returned in *DIR_UPDATE_PTR (used for updating predictor state),
   and the non-speculative top-of-stack is returned in stack_recover_idx 
   (used for recovering ret-addr stack after mis-predict).  */
md_addr_t				/* predicted branch target addr */
bpred_lookup(struct bpred_t *pred,	/* branch predictor instance */
	     md_addr_t baddr,		/* branch address */
	     md_addr_t btarget,		/* branch target if taken */
	     enum md_opcode op,		/* opcode of instruction */
	     int is_call,		/* non-zero if inst is fn call */
	     int is_return,		/* non-zero if inst is fn return */
	     struct bpred_update_t *dir_update_ptr, /* pred state pointer */
	     int *stack_recover_idx,	/* Non-speculative top-of-stack; used on mispredict recovery */
		 int flow_mode)  /* Flow mode (0=FWD; 1=REV) */
{
	struct bpred_btb_ent_t *pbtb = NULL;
	int index, i;
	bool_t invert = FALSE;
	  
	int fwd_valid_outcome = NULL;
	int rev_valid_outcome = NULL;

	if (!dir_update_ptr)
		panic("no bpred update record");

	/* if this is not a branch, return not-taken */
	if (!(MD_OP_FLAGS(op) & F_CTRL))
		return 0;
   
   // REV addresses; 
   // These are supposed to be opposites of the FWD but 
   // without ISA updates and knowledge of addrs and targets, 
   // these need to be the same.
   md_addr_t rbaddr = baddr;		/* REV Branch ADDR; for now, setting to baddr */
   md_addr_t rbtarget = btarget;	/* REV Branch Target; for now, setting to btarget */
   

	if (!flow_mode) {
		pred->lookups++;
	} else {
		pred->reverse_lookups++;
	}

	// FWD DIR Pointers
	dir_update_ptr->fwd_dir.ras = FALSE;
	dir_update_ptr->fwd_pdir1 = NULL;
	dir_update_ptr->fwd_pdir2 = NULL;
	dir_update_ptr->fwd_pmeta = NULL;
	
	// REV DIR Pointers
	dir_update_ptr->rev_dir.ras = FALSE;
	dir_update_ptr->rev_pdir1 = NULL;
	dir_update_ptr->rev_pdir2 = NULL;
	dir_update_ptr->rev_pmeta = NULL;
	
	/* Except for jumps, get a pointer to direction-prediction bits */
	switch (pred->class) {
		case BPredComb:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				char *bimod, *twolev, *meta;
			  
				bimod = bpred_dir_lookup (pred->fwd_dirpred.bimod, baddr);
				twolev = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);
				meta = bpred_dir_lookup (pred->fwd_dirpred.meta, baddr);
				  
				dir_update_ptr->fwd_pmeta = meta;
				dir_update_ptr->fwd_dir.meta  = (*meta >= 2);
				dir_update_ptr->fwd_dir.bimod = (*bimod >= 2);
				dir_update_ptr->fwd_dir.twolev  = (*twolev >= 2);
				  
				if (*meta >= 2) {
					dir_update_ptr->fwd_pdir1 = twolev;
					dir_update_ptr->fwd_pdir2 = bimod;
				} else {
					dir_update_ptr->fwd_pdir1 = bimod;
					dir_update_ptr->fwd_pdir2 = twolev;
				}
			}
			break;
			
		case BPred2Level:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if FHB has valid past outcome results
				if (pred->fhb.fv[pred->fhb.top]) {
					fwd_valid_outcome = pred->fhb.o[pred->fhb.top];
				}
					
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);
			}
			break;
			
		case BPredTSBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->fwd_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
				
				/* incr head but prevent from going out of bounds*/
				if (pred->fwd_dirpred.tsbp->ts.head >= pred->fwd_dirpred.tsbp->ts.correctness_width) {
					pred->fwd_dirpred.tsbp->ts.head = 0;
				} else {
					pred->fwd_dirpred.tsbp->ts.head++;
				}

				/*if in replay mode and corretness buffer head indicates base predictor mistake*/
				if(pred->fwd_dirpred.tsbp->ts.replay 
					&& pred->fwd_dirpred.tsbp->ts.enabled 
					&& (pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.head] == 0)) {
					invert = TRUE; 
				}
			}
			break;
			
		case BPredCHBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->fwd_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size

				/*if enabled, replay bit set, correctness bits is 0, and src_pc matches baddr predictor is inverted*/
				if (pred->fwd_dirpred.chbp->chbp.enabled 
					&& pred->fwd_dirpred.chbp->chbp.cht_replay[key] 				// Only perform correction if replay is on
					&& !pred->fwd_dirpred.chbp->chbp.cht_correct[key] 			// Check for past correctness history
					&& (pred->fwd_dirpred.chbp->chbp.cht_spc[key] == baddr)) { 	// Check stored source pc is same as baddr
					invert = TRUE; 
				}
			}
			break;
			
		case BPredOB:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if OB or FHB have valid past outcome results
				if (pred->ob.fv[pred->ob.end]) {
					fwd_valid_outcome = pred->ob.oc[pred->ob.end];
				} else if (pred->fhb.fv[pred->fhb.top]) {
					fwd_valid_outcome = pred->fhb.o[pred->fhb.top];
				}
					
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);
			}
			break;
			
		case BPredOHT:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if FHB or OHT have valid past outcome results
				if (pred->fhb.fv[pred->fhb.top]) {
					fwd_valid_outcome = pred->fhb.o[pred->fhb.top];
				} else {
					int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
					key = key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on predictor table size
					
					if (pred->fwd_dirpred.oht->oht.valid[key]) {
						fwd_valid_outcome = pred->fwd_dirpred.oht->oht.oc[key];
						
						// Only valid once in the OHT when in FWD mode
						if (!flow_mode) {
							pred->fwd_dirpred.oht->oht.valid[key] = 0;
						}
					}
				}
					
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);
			}
			break;
			
		case BPredMBP:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if OB or OHT have valid past outcome results
				if (pred->ob.fv[pred->ob.end]) {
					fwd_valid_outcome = pred->ob.oc[pred->ob.end];
				} else {
					int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
					key = key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on predictor table size
					
					if (pred->fwd_dirpred.oht->oht.valid[key]) {
						fwd_valid_outcome = pred->fwd_dirpred.oht->oht.oc[key];
						
						// Only valid once in the OHT when in FWD mode
						if (!flow_mode) {
							pred->fwd_dirpred.oht->oht.valid[key] = 0;
						}
					}
				}
					
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup (pred->fwd_dirpred.twolev, baddr);
			}
			break;
			
		case BPred2bit:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->fwd_pdir1 =
				bpred_dir_lookup (pred->fwd_dirpred.bimod, baddr);
			}
			break;
			
		case BPredTaken:
			return btarget;
			
		case BPredNotTaken:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				return baddr + sizeof(md_inst_t);
			} else {
				return btarget;
			}
			
		default:
			panic("bogus predictor class");
	}
	
	// REV Predictors
	  
	/* Except for jumps, get a pointer to direction-prediction bits */
	switch (pred->class) {
		case BPredComb:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				char *bimod, *twolev, *meta;
				
				bimod = bpred_dir_lookup (pred->rev_dirpred.bimod, baddr);
				twolev = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);
				meta = bpred_dir_lookup (pred->rev_dirpred.meta, baddr);
				
				dir_update_ptr->rev_pmeta = meta;
				dir_update_ptr->rev_dir.meta  = (*meta >= 2);
				dir_update_ptr->rev_dir.bimod = (*bimod >= 2);
				dir_update_ptr->rev_dir.twolev  = (*twolev >= 2);
				
				if (*meta >= 2) {
					dir_update_ptr->rev_pdir1 = twolev;
					dir_update_ptr->rev_pdir2 = bimod;
				} else {
					dir_update_ptr->rev_pdir1 = bimod;
					dir_update_ptr->rev_pdir2 = twolev;
				}
			}
			break;
			
		case BPred2Level:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if FHB has valid past outcome results
				if (pred->fhb.rv[pred->fhb.bot]) {
					rev_valid_outcome = pred->fhb.o[pred->fhb.bot];
				}
					
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);
			}
			break;
			
		case BPredTSBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->rev_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->rev_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
				
				/* incr head but prevent from going out of bounds*/
				if (pred->rev_dirpred.tsbp->ts.head >= pred->rev_dirpred.tsbp->ts.correctness_width) {
					pred->rev_dirpred.tsbp->ts.head = 0;
				} else {
					pred->rev_dirpred.tsbp->ts.head++;
				}

				/*if in replay mode and corretness buffer head indicates base predictor mistake*/
				if(pred->rev_dirpred.tsbp->ts.replay 
					&& pred->rev_dirpred.tsbp->ts.enabled 
					&& (pred->rev_dirpred.tsbp->ts.correctness_buffer[pred->rev_dirpred.tsbp->ts.head] == 0)) {
					invert = TRUE; 
				}
			}
			break;
			
		case BPredCHBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->rev_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->rev_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size

				/*if enabled, replay bit set, correctness bits is 0, and src_pc matches baddr predictor is inverted*/
				if (pred->rev_dirpred.chbp->chbp.enabled 
					&& pred->rev_dirpred.chbp->chbp.cht_replay[key] 				// Only perform correction if replay is on
					&& !pred->rev_dirpred.chbp->chbp.cht_correct[key] 			// Check for past correctness history
					&& (pred->rev_dirpred.chbp->chbp.cht_spc[key] == baddr)) { 	// Check stored source pc is same as baddr
					invert = TRUE; 
				}
			}
			break;
			
		case BPredOB:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if OB or FHB have valid past outcome results
				if (pred->ob.rv[pred->ob.beg]) {
					rev_valid_outcome = pred->ob.oc[pred->ob.beg];
				} else if (pred->fhb.rv[pred->fhb.bot]) {
					rev_valid_outcome = pred->fhb.o[pred->fhb.bot];
				}
					
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);
			}
			break;
			
		case BPredOHT:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if FHB or OHT have valid past outcome results
				if (pred->fhb.rv[pred->fhb.bot]) {
					rev_valid_outcome = pred->fhb.o[pred->fhb.bot];
				} else {
					int key = key_from_features (pred->rev_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
					key = key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on predictor table size
					
					if (pred->rev_dirpred.oht->oht.valid[key]) {
						rev_valid_outcome = pred->rev_dirpred.oht->oht.oc[key];
						
						// Only valid once in the OHT when in REV mode
						if (flow_mode) {
							pred->rev_dirpred.oht->oht.valid[key] = 0;
						}
					}
				}
					
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);
			}
			break;
			
		case BPredMBP:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if OB or OHT have valid past outcome results
				if (pred->ob.rv[pred->ob.beg]) {
					rev_valid_outcome = pred->ob.oc[pred->ob.beg];
				} else {
					int key = key_from_features (pred->rev_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
				
					key = key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on predictor table size
					
					if (pred->rev_dirpred.oht->oht.valid[key]) {
						rev_valid_outcome = pred->rev_dirpred.oht->oht.oc[key];
						
						// Only valid once in the OHT when in REV mode
						if (flow_mode) {
							pred->rev_dirpred.oht->oht.valid[key] = 0;
						}
					}
				}
					
				dir_update_ptr->rev_pdir1 = bpred_dir_lookup (pred->rev_dirpred.twolev, baddr);
			}
			break;
			
		case BPred2bit:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->rev_pdir1 =
				bpred_dir_lookup (pred->rev_dirpred.bimod, baddr);
			}
			break;
			
		case BPredTaken:
			return btarget;
			
		case BPredNotTaken:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				return baddr + sizeof(md_inst_t);
			} else {
				return btarget;
			}
			
		default:
			panic("bogus predictor class");
	}
	
	
   /*
	* We have a stateful predictor, and have gotten a pointer into the
	* direction predictor (except for jumps, for which the ptr is null)
	*
	* record pre-pop TOS; if this branch is executed speculatively
	* and is squashed, we'll restore the TOS and hope the data
	* wasn't corrupted in the meantime.
	*
	* For REV mode, return requires push and call requires pop 
	*
	*/
	
	if (!flow_mode) {
		if (pred->retstack.size)
			*stack_recover_idx = pred->retstack.tos;
		else
			*stack_recover_idx = 0;

		/* if this is a return, pop return-address stack */
		if (is_return && pred->retstack.size) {
			md_addr_t target = pred->retstack.stack[pred->retstack.tos].target;
			
			pred->retstack.tos = (pred->retstack.tos + pred->retstack.size - 1) % pred->retstack.size;
			pred->retstack_pops++;
			
			dir_update_ptr->fwd_dir.ras = TRUE; /* using RAS here */
			
			return target;
		}

#ifndef RAS_BUG_COMPATIBLE
		/* if function call, push return-address onto return-address stack */
		if (is_call && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = baddr + sizeof(md_inst_t);
			pred->retstack_pushes++;
		}
#endif /* !RAS_BUG_COMPATIBLE */
	} else {
		if (pred->retstack.size)
			*stack_recover_idx = pred->retstack.tos;
		else
			*stack_recover_idx = 0;

		/* Since this is reversed, return constitutes a push and call a pop  */
		if (is_return && pred->retstack.size) {
#ifndef RAS_BUG_COMPATIBLE
			/* if this is a return, push btarget onto return-address stack */
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = btarget;
			
			pred->reverse_retstack_pushes++;
#endif /* !RAS_BUG_COMPATIBLE */
			return btarget;
		}
		
		/* if function call, pop return-address stack */
		if (is_call && pred->retstack.size) {
			// Don't need target for this as there is no reversible ISA yet implemented
			//md_addr_t target = pred->retstack.stack[pred->retstack.tos].target;
			
			pred->retstack.tos = (pred->retstack.tos + pred->retstack.size - 1) % pred->retstack.size;
			pred->reverse_retstack_pops++;
			
			dir_update_ptr->rev_dir.ras = TRUE; /* using RAS here */
			
			//return target;
		}
	}
  
	/* Get a pointer into the BTB early as we may need target and it's already set after FWD mode */
	index = (baddr >> MD_BR_SHIFT) & (pred->btb.sets - 1);

	if (pred->btb.assoc > 1) {
		index *= pred->btb.assoc;

		/* Now we know the set; look for a PC match */
		for (i = index; i < (index+pred->btb.assoc) ; i++)
			if (pred->btb.btb_data[i].addr == baddr) {
				/* match */
				pbtb = &pred->btb.btb_data[i];
				break;
			}
	} else {
		pbtb = &pred->btb.btb_data[index];
		
		if (pbtb->addr != baddr)
			pbtb = NULL;
	}

  /*
   * We now also have a pointer into the BTB for a hit, or NULL otherwise
   */

	/* if this is a jump, ignore predicted direction; we know it's taken. */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) == (F_CTRL|F_UNCOND)) {
		return (pbtb ? pbtb->target : 1);
	}
	
	int prediction;
		
	/* otherwise we have a conditional branch */
	if (!flow_mode) {
		// There was a fwd valid outcome
		if (fwd_valid_outcome != NULL) {
			prediction = fwd_valid_outcome;
		} else {
			prediction = (*(dir_update_ptr->fwd_pdir1) >= 2);
			
			if (invert) {
				prediction = !prediction;
			}
		}
	} else {
		// There was a rev valid outcome
		if (rev_valid_outcome != NULL) {
			prediction = rev_valid_outcome;
		} else {
			prediction = (*(dir_update_ptr->rev_pdir1) >= 2);
			
			if (invert) {
				prediction = !prediction;
			}
		}
	}
	  
	md_addr_t ret_target; 
	
	if (pbtb == NULL) {
		/* BTB miss -- just return a predicted direction */
		ret_target = 1;
	} else {
		/* BTB hit, so return target if it's a predicted-taken branch */
		ret_target = pbtb->target;
	}
	
	return (prediction
		? /* taken */ ret_target
		: /* not taken */ 0
	);
}

/* Speculative execution can corrupt the ret-addr stack.  So for each
 * lookup we return the top-of-stack (TOS) at that point; a mispredicted
 * branch, as part of its recovery, restores the TOS using this value --
 * hopefully this uncorrupts the stack. */
void
bpred_recover(struct bpred_t *pred,	/* branch predictor instance */
	      md_addr_t baddr,		/* branch address */
	      int stack_recover_idx)	/* Non-speculative top-of-stack;
					 * used on mispredict recovery */
{
  if (pred == NULL)
    return;

  pred->retstack.tos = stack_recover_idx;
}

/* update the branch predictor, only useful for stateful predictors; updates
   entry for instruction type OP at address BADDR.  BTB only gets updated
   for branches which are taken.  Inst was determined to jump to
   address BTARGET and was taken if TAKEN is non-zero.  Predictor 
   statistics are updated with result of prediction, indicated by CORRECT and 
   PRED_TAKEN, predictor state to be updated is indicated by *DIR_UPDATE_PTR 
   (may be NULL for jumps, which shouldn't modify state bits).  Note if
   bpred_update is done speculatively, branch-prediction may get polluted. */
void
bpred_update(struct bpred_t *pred,	/* branch predictor instance */
	     md_addr_t baddr,		/* branch address */
	     md_addr_t btarget,		/* resolved branch target */
		 int taken,			/* non-zero if branch was taken */
	     int pred_taken,		/* non-zero if branch was pred taken */
	     int correct,		/* was earlier addr prediction ok? */
	     enum md_opcode op,		/* opcode of instruction */
	     struct bpred_update_t *dir_update_ptr,  /* FWD pred state pointer */
		 int flow_mode)  /* Flow mode (0=FWD; 1=REV) */
{
	struct bpred_btb_ent_t *pbtb = NULL;
	struct bpred_btb_ent_t *lruhead = NULL, *lruitem = NULL;
	int index, i;
  
   /* don't change bpred state for non-branch instructions or if this
	* is a stateless predictor*/
	if (!(MD_OP_FLAGS(op) & F_CTRL))
		return;

	/* Have a branch here */

	/* Now accounting for FWD or REV stats */
	if (!flow_mode) {
		if (correct)
			pred->addr_hits++;

		if (!!pred_taken == !!taken)
			pred->dir_hits++;
		else
			pred->misses++;

		if (dir_update_ptr->fwd_dir.ras) {
			pred->used_ras++;
			
			if (correct)
				pred->ras_hits++;
		} else if ((MD_OP_FLAGS(op) & (F_CTRL|F_COND)) == (F_CTRL|F_COND)) {
			if (dir_update_ptr->fwd_dir.meta)
				pred->used_2lev++;
			else
				pred->used_bimod++;
		}
		
		/* keep stats about JR's; also, but don't change any bpred state for JR's
		* which are returns unless there's no retstack */
		if (MD_IS_INDIR(op)) {
			pred->jr_seen++;
			if (correct)
				pred->jr_hits++;
		  
			if (!dir_update_ptr->fwd_dir.ras) {
				pred->jr_non_ras_seen++;
				
				if (correct)
					pred->jr_non_ras_hits++;
			} else {
				/* return that used the ret-addr stack; no further work to do */
				return;
			}
		}
	} else {
		if (correct)
			pred->reverse_addr_hits++;

		if (!!pred_taken == !!taken)
			pred->reverse_dir_hits++;
		else
			pred->reverse_misses++;

		if (dir_update_ptr->rev_dir.ras) {
			pred->reverse_used_ras++;
			
			if (correct)
				pred->reverse_ras_hits++;
		} else if ((MD_OP_FLAGS(op) & (F_CTRL|F_COND)) == (F_CTRL|F_COND)) {
			if (dir_update_ptr->rev_dir.meta)
				pred->reverse_used_2lev++;
			else
				pred->reverse_used_bimod++;
		}
		
		/* keep stats about rev JR's (CALLs); also, but don't change any bpred state for JR's
		* which are returns unless there's no retstack */
		if (MD_IS_CALL(op)) {
			pred->reverse_jr_seen++;
			
			if (correct)
				pred->reverse_jr_hits++;
		  
			if (!dir_update_ptr->rev_dir.ras) {
				pred->reverse_jr_non_ras_seen++;
				
				if (correct)
					pred->reverse_jr_non_ras_hits++;
			} else {
				/* return that used the ret-addr stack; no further work to do */
				return;
			}
		}
	}

	/* Can exit now if this is a stateless predictor */
	if (pred->class == BPredNotTaken || pred->class == BPredTaken)
		return;

   /* 
	* Now we know the branch didn't use the ret-addr stack, and that this
	* is a stateful predictor 
	*
	* Added handling for REV where return constitutes a push
	*/
 
#ifdef RAS_BUG_COMPATIBLE
	if (!flow_mode) {
		/* if function call, push return-address onto return-address stack */
		if (MD_IS_CALL(op) && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = baddr + sizeof(md_inst_t);
		
			pred->retstack_pushes++;
		}
	} else {
		if (MD_IS_INDIR(op) && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = btarget;
		
			pred->reverse_retstack_pushes++;
		}
	}
#endif /* RAS_BUG_COMPATIBLE */

	// Updating BTB earlier to get target for rev pred updates
	
	/* find BTB entry if it's a taken branch (don't allocate for non-taken) */
	if (taken) {
		index = (baddr >> MD_BR_SHIFT) & (pred->btb.sets - 1);
      
		if (pred->btb.assoc > 1) {
			index *= pred->btb.assoc;
	  
		   /* Now we know the set; look for a PC match; also identify
			* MRU and LRU items */
			for (i = index; i < (index+pred->btb.assoc) ; i++) {
				if (pred->btb.btb_data[i].addr == baddr) {
					/* match */
					assert(!pbtb);
					pbtb = &pred->btb.btb_data[i];
				}
	      
				dassert(pred->btb.btb_data[i].prev != pred->btb.btb_data[i].next);
				if (pred->btb.btb_data[i].prev == NULL) {
					/* this is the head of the lru list, ie current MRU item */
					dassert(lruhead == NULL);
					lruhead = &pred->btb.btb_data[i];
				}
				
				if (pred->btb.btb_data[i].next == NULL) {
					/* this is the tail of the lru list, ie the LRU item */
					dassert(lruitem == NULL);
					lruitem = &pred->btb.btb_data[i];
				}
			}
			dassert(lruhead && lruitem);
	  
			if (!pbtb)
				/* missed in BTB; choose the LRU item in this set as the victim */
				pbtb = lruitem;	
				/* else hit, and pbtb points to matching BTB entry */
			  
		   /* Update LRU state: selected item, whether selected because it
			* matched or because it was LRU and selected as a victim, becomes 
			* MRU */
			if (pbtb != lruhead) {
				/* this splices out the matched entry... */
				if (pbtb->prev)
					pbtb->prev->next = pbtb->next;
				if (pbtb->next)
					pbtb->next->prev = pbtb->prev;
				/* ...and this puts the matched entry at the head of the list */
				pbtb->next = lruhead;
				pbtb->prev = NULL;
				lruhead->prev = pbtb;
				dassert(pbtb->prev || pbtb->next);
				dassert(pbtb->prev != pbtb->next);
			}
			/* else pbtb is already MRU item; do nothing */
		} else {
			pbtb = &pred->btb.btb_data[index];
		}
      
	   /* 
		* Now 'p' is a possibly null pointer into the direction prediction table, 
		* and 'pbtb' is a possibly null pointer into the BTB (either to a 
		* matched-on entry or a victim which was LRU in its set)
		*/
		   
		/* update BTB (but only for taken branches) */
		if (pbtb) {
			/* update current information */
			dassert(taken);

			if (pbtb->addr == baddr) {
				if (!correct)
					pbtb->target = btarget;
			} else {
				/* enter a new branch in the table */
				pbtb->addr = baddr;
				pbtb->op = op;
				pbtb->target = btarget;
			}
		}
	}

	int fwd_key, rev_key, fv, rv, oc;
	md_addr_t fhb_addr;
  

	/* Get keys before updating L1 table*/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && ((pred->class == BPredTSBP) || (pred->class == BPredCHBP) || (pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		// Get keys before updating GHR!!!!
		fwd_key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
		rev_key = key_from_features (pred->rev_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
	}

	/* update L1 table if appropriate */
	/* L1 table is updated unconditionally for combining predictor too */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) &&
	((pred->class == BPredComb) || (pred->class == BPred2Level) || (pred->class == BPredTSBP) || (pred->class == BPredCHBP) || (pred->class == BPredOB) || (pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		int l1index, shift_reg;
      
		/* also update appropriate L1 history registers; Shift R or L depending on flow mode */
		if (!flow_mode) {
			l1index = (baddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
			shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
			pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);
			
			l1index = (baddr >> MD_BR_SHIFT) & (pred->rev_dirpred.twolev->config.two.l1size - 1);
			shift_reg = (pred->rev_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
			pred->rev_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->rev_dirpred.twolev->config.two.shift_width) - 1);
		} else {
			l1index = (baddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
			shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] >> 1) | ((!!taken) << (pred->fwd_dirpred.twolev->config.two.shift_width - 1));
			pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);
			
			l1index = (baddr >> MD_BR_SHIFT) & (pred->rev_dirpred.twolev->config.two.l1size - 1);
			shift_reg = (pred->rev_dirpred.twolev->config.two.shiftregs[l1index] >> 1) | ((!!taken) << (pred->rev_dirpred.twolev->config.two.shift_width - 1));
			pred->rev_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->rev_dirpred.twolev->config.two.shift_width) - 1);
		}
	
		//Update FHB and shift up or down depending on flow mode
		if (!flow_mode) {
			fhb_addr = pred->fhb.addr[pred->fhb.top];
			fv = pred->fhb.fv[pred->fhb.top];
			rv = pred->fhb.rv[pred->fhb.top];
			oc = pred->fhb.o[pred->fhb.top];
			
			if (pred->fhb.top > 0) {
				pred->fhb.top--;
			} else {
				pred->fhb.top = pred->fhb.size - 1;
			}
			
			if (pred->fhb.bot > 0) {
				pred->fhb.bot--;
			} else {
				pred->fhb.bot = pred->fhb.size - 1;
			}
			
			if (pbtb) {
				pred->fhb.addr[pred->fhb.bot] = pbtb->target;
			} else {
				pred->fhb.addr[pred->fhb.bot] = NULL;
			}
			
			pred->fhb.fv[pred->fhb.bot] = 0;
			pred->fhb.rv[pred->fhb.bot] = 1;
			pred->fhb.o[pred->fhb.bot] = taken;
		} else {
			fhb_addr = pred->fhb.addr[pred->fhb.bot];
			fv = pred->fhb.fv[pred->fhb.bot];
			rv = pred->fhb.rv[pred->fhb.bot];
			oc = pred->fhb.o[pred->fhb.bot];
			
			if (pred->fhb.top < (pred->fhb.size - 1)) {
				pred->fhb.top++;
			} else {
				pred->fhb.top = 0;
			}
			
			if (pred->fhb.bot < (pred->fhb.size - 1)) {
				pred->fhb.bot++;
			} else {
				pred->fhb.bot = 0;
			}
			
			if (pbtb) {
				pred->fhb.addr[pred->fhb.top] = pbtb->addr;
			} else {
				pred->fhb.addr[pred->fhb.top] = NULL;
			}
			
			pred->fhb.fv[pred->fhb.top] = 1;
			pred->fhb.rv[pred->fhb.top] = 0;
			pred->fhb.o[pred->fhb.top] = taken;
		}
    }

	/***********************IF TS, also update correctness buffer*****************************/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredTSBP)) {
		//int l1index, shift_reg;
		
		// Get key before updating GHR!!!!
		//int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
    
		/* update L1 table, same as 2lev/comb predictors above this */
		//l1index = (baddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
		//shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
		//pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);

		int base_outcome;
		int ts_outcome;
		//unsigned int key;  /*added typedef in tsbp.h file*/
      	
		fwd_key = fwd_key & (pred->fwd_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
		
		/*Set the base outcome; Also if in replay mode and ts_outcome is incorrect, turn off replay mode*/
		if (pred->fwd_dirpred.tsbp->ts.replay) {
			pred->replays++;
			ts_outcome = pred_taken;

			if (pred->fwd_dirpred.tsbp->ts.enabled && (pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.head] == 0)) {
				base_outcome = !pred_taken;
			} else {
				base_outcome = pred_taken;
			}

			if (!!ts_outcome != !!taken) {
				pred->fwd_dirpred.tsbp->ts.replay = FALSE;
			}
		} else {
			base_outcome = pred_taken;
		}
		
		/* incr tail but prevent from going out of bounds*/
		if (pred->fwd_dirpred.tsbp->ts.tail >= pred->fwd_dirpred.tsbp->ts.correctness_width) {
			pred->fwd_dirpred.tsbp->ts.tail = 0;
		} else {
			pred->fwd_dirpred.tsbp->ts.tail++;
		}
		
		/*determine if actual outcome (taken) of predicted direction is correct and update correctness buffer*/
		pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.tail] = (!!base_outcome == !!taken);  /*1 = base predictor correct. 0 = prediction incorrect*/
      
		/*if incorrect base prediction, update head table*/
		if (!!base_outcome != !!taken) {
			if(!pred->fwd_dirpred.tsbp->ts.replay && pred->fwd_dirpred.tsbp->ts.head_table[fwd_key] != NULL) { /*if not in replay mode, update head and set replay flag*/
				pred->fwd_dirpred.tsbp->ts.head = pred->fwd_dirpred.tsbp->ts.head_table[fwd_key];
				pred->fwd_dirpred.tsbp->ts.replay = TRUE;
			}
		}

		pred->fwd_dirpred.tsbp->ts.head_table[fwd_key] = pred->fwd_dirpred.tsbp->ts.tail;   //else update head table
	}
  
	/***********************IF CHBP, also update correctness history table*****************************/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredCHBP)) {
		//int l1index, shift_reg;
		
		// Get key before updating GHR!!!!
		//int key = key_from_features (pred->fwd_dirpred.twolev, baddr); // Get unmasked key from GHR and PC
		
		/* update L1 table, same as 2lev/comb predictors above this */
		//l1index = (baddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
		//shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
		//pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);

		int base_outcome;
		int chbp_outcome;
		
		fwd_key = fwd_key & (pred->fwd_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size
		
		bool_t just_disabled_replay = FALSE;

		/*Set the base outcome; Also if in replay mode and chbp_outcome is incorrect, turn off replay mode*/
		if (pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] && (pred->fwd_dirpred.chbp->chbp.cht_spc[fwd_key] == baddr)) {
			pred->replays++;
			chbp_outcome = pred_taken;

			if (pred->fwd_dirpred.chbp->chbp.enabled && !pred->fwd_dirpred.chbp->chbp.cht_correct[fwd_key]) {
				base_outcome = !pred_taken;
			} else {
				base_outcome = pred_taken;
			}

			if (!!chbp_outcome != !!taken) {
				pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] = FALSE;
				just_disabled_replay = TRUE;
			}
		} else {
			base_outcome = pred_taken;
		}
		
		/*determine if actual outcome (taken) of predicted direction is correct and update correctness bit*/
		pred->fwd_dirpred.chbp->chbp.cht_correct[fwd_key] = (!!base_outcome == !!taken);  /*1 = base predictor correct. 0 = prediction incorrect*/
		  
		/*if not in replay mode and base outcome incorrect, set replay flag as long as it wasn't just disabled*/
		if ((!!base_outcome != !!taken) && !pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] && !just_disabled_replay) {
			pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] = TRUE;
		}
		
		/* Update other correctness history table bits */
		pred->fwd_dirpred.chbp->chbp.cht_dpc[fwd_key] = btarget;
		pred->fwd_dirpred.chbp->chbp.cht_valid[fwd_key] = TRUE;
		pred->fwd_dirpred.chbp->chbp.cht_spc[fwd_key] = baddr;
		
		/* make sure no other matching spc keys are valid*/
		//int valid_key;
		
		//for (valid_key = 0; valid_key < pred->fwd_dirpred.chbp->chbp.cht_size; valid_key++) {
		//	if ((pred->fwd_dirpred.chbp->chbp.cht_spc[valid_key] == baddr) && pred->fwd_dirpred.chbp->chbp.cht_valid[valid_key] && (valid_key != key)) {
		//		pred->fwd_dirpred.chbp->chbp.cht_valid[valid_key] = FALSE;
		//	}
		//}
		
	}
	
	/* Update OB if OB or MBP */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && ((pred->class == BPredOB) || (pred->class == BPredMBP))) {
		//Update OB and shift up (right) or down (left) depending on flow mode
		if (!flow_mode) {
			if (pred->ob.end > 0) {
				pred->ob.end--;
			} else {
				pred->ob.end = pred->ob.width - 1;
			}
			
			if (pred->ob.beg > 0) {
				pred->ob.beg--;
			} else {
				pred->ob.beg = pred->ob.width - 1;
			}
			
			pred->ob.fv[pred->ob.beg] = 0;
			pred->ob.rv[pred->ob.beg] = 1;
			pred->ob.oc[pred->ob.beg] = taken;
		} else {
			if (pred->ob.end < (pred->ob.width - 1)) {
				pred->ob.end++;
			} else {
				pred->ob.end = 0;
			}
			
			if (pred->ob.beg < (pred->ob.width - 1)) {
				pred->ob.beg++;
			} else {
				pred->ob.beg = 0;
			}
			
			pred->ob.fv[pred->ob.end] = 1;
			pred->ob.rv[pred->ob.end] = 0;
			pred->ob.oc[pred->ob.end] = taken;
		}
	}
	
	/* Update OHT if OHT or MBP */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && ((pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		// Update only one OHT depending on the flow and only if fhb addr is set
		if (fhb_addr) {
			fwd_key = key_from_features (pred->fwd_dirpred.twolev, fhb_addr); // Get unmasked key from GHR and PC
			rev_key = key_from_features (pred->rev_dirpred.twolev, fhb_addr); // Get unmasked key from GHR and PC
			
			fwd_key = fwd_key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
			rev_key = rev_key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
			
			if (!flow_mode) {
				pred->rev_dirpred.oht->oht.oc[rev_key] = oc; //Update with oc since it is the past outcome that was waiting in the FHB
				pred->rev_dirpred.oht->oht.valid[rev_key] = 1;
			} else {
				pred->fwd_dirpred.oht->oht.oc[fwd_key] = oc;
				pred->fwd_dirpred.oht->oht.valid[fwd_key] = 1;
			}
		}		
	}

	/* update state (but not for jumps) */
	if (dir_update_ptr->fwd_pdir1) {  
        if (taken) {
			if (*dir_update_ptr->fwd_pdir1 < 3)
				++*dir_update_ptr->fwd_pdir1;
		} else { /* not taken */
			if (*dir_update_ptr->fwd_pdir1 > 0)
				--*dir_update_ptr->fwd_pdir1;
        }
    }
	
	/* Acounting for reverse dir pointer now */
	if (dir_update_ptr->rev_pdir1) {  
        if (taken) {
			if (*dir_update_ptr->rev_pdir1 < 3)
				++*dir_update_ptr->rev_pdir1;
		} else { /* not taken */
			if (*dir_update_ptr->rev_pdir1 > 0)
				--*dir_update_ptr->rev_pdir1;
        }
    }

	/* combining predictor also updates second predictor and meta predictor */
	/* second direction predictor */
	if (dir_update_ptr->fwd_pdir2) {
		if (taken) {
			if (*dir_update_ptr->fwd_pdir2 < 3)
				++*dir_update_ptr->fwd_pdir2;
		} else { /* not taken */
			if (*dir_update_ptr->fwd_pdir2 > 0)
				--*dir_update_ptr->fwd_pdir2;
		}
    }
	
	/* Acounting for reverse dir pointer now */
	if (dir_update_ptr->rev_pdir2) {  
        if (taken) {
			if (*dir_update_ptr->rev_pdir2 < 3)
				++*dir_update_ptr->rev_pdir2;
		} else { /* not taken */
			if (*dir_update_ptr->rev_pdir2 > 0)
				--*dir_update_ptr->rev_pdir2;
        }
    }

	/* meta predictor */
	if (dir_update_ptr->fwd_pmeta) {
		if (dir_update_ptr->fwd_dir.bimod != dir_update_ptr->fwd_dir.twolev) {
			/* we only update meta predictor if directions were different */
			if (dir_update_ptr->fwd_dir.twolev == (unsigned int)taken) {
				/* 2-level predictor was correct */
				if (*dir_update_ptr->fwd_pmeta < 3)
					++*dir_update_ptr->fwd_pmeta;
			} else {
				/* bimodal predictor was correct */
				if (*dir_update_ptr->fwd_pmeta > 0)
					--*dir_update_ptr->fwd_pmeta;
			}
		}
    }
	
	/* Acounting for reverse dir pointer now */
	if (dir_update_ptr->rev_pmeta) {
		if (dir_update_ptr->rev_dir.bimod != dir_update_ptr->rev_dir.twolev) {
			/* we only update meta predictor if directions were different */
			if (dir_update_ptr->rev_dir.twolev == (unsigned int)taken) {
				/* 2-level predictor was correct */
				if (*dir_update_ptr->rev_pmeta < 3)
					++*dir_update_ptr->rev_pmeta;
			} else {
				/* bimodal predictor was correct */
				if (*dir_update_ptr->rev_pmeta > 0)
					--*dir_update_ptr->rev_pmeta;
			}
		}
    }
}
