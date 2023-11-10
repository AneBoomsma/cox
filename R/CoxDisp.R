
NumToString <- function(x, dec) {

  #Return a number as string with n decimal places
  string <- formatC(x, digits = dec, format = 'f')
  return(string)

}

PVal <- function (x) {

  #Return <0.001 if p is less than 0.001
  if(is.na(x)){
    return("NA")
  } else {
    if (x < 0.001) {
      return('<0.001')
    } else {
      return(NumToString(x, 3))
    }
  }

}

.formatCoef <- function(exp.coef, coef.lower, coef.upper) {

  #Return formatted coefficients with 95% CI
  coef.string <- paste0(NumToString(exp.coef, 2),
                        ' [',
                        NumToString(coef.lower, 2),
                        ', ',
                        NumToString(coef.upper, 2),
                        ']')
  return(coef.string)

}

CoxSummary <- function(models = models, var = var, model.names = NULL, summ = TRUE, styled = TRUE) {

  output.df <- data.frame()

  #Which model in the loop
  model.n <- 0

  for (model in models) {

    independent.var <- names(model$coefficients)[match(TRUE, grepl(var,names(model$coefficients)))]

    #Calculate model number it is
    model.n <- model.n + 1

    #Summary of the Cox Model
    model.summ <- summary(model)

    #Get variable names in the model
    variables.vector <- as.vector(names(model$coefficients))

    #Get variable names from the cox.zph function - this can differ from
    #variable names in the model in case of categorical variables
    vector.var.in.cox.zph <- row.names(cox.zph(model)$table)
    for (var.name in vector.var.in.cox.zph) {
      schoenfeld.res.test.p <- cox.zph(model)$table[row.names(cox.zph(model)$table)[match(TRUE, grepl(var.name, row.names(cox.zph(model)$table)))],'p']
      #Check if the model meets the proportional hazards assumption
      if (schoenfeld.res.test.p < 0.05) {
        warning(paste0('Proportional hazard assumption might not be reasonable',
                     '\n\nThe variable:',
                     var.name,
                     ' in model with the variables:\n',
                     paste0(vector.var.in.cox.zph[-length(vector.var.in.cox.zph)], collapse = ', '),
                     '\ndoes not meet the proportional hazards assumption.',
                     '\n\nSchoenfeld residuals based statistical test: p = ',
                     round(schoenfeld.res.test.p, 2),
                     ' \n\n'))
        }
      }

    #Model is named as model + n
    model.name <- paste0('Model ', model.n)

    #If summ is TRUE, only a summary of coefficients of the independent var will be returned
    if (summ == FALSE) {

      var.n <- 0

      for (var in variables.vector) {

        var.n <- var.n + 1

        #P value
        p.val <- PVal(model.summ$coefficients[var, 'Pr(>|z|)'])

        #Get Exp(coef) and 95% CI from the model
        exp.coef <- model.summ$conf.int[var, 'exp(coef)']
        coef.lower <- model.summ$conf.int[var, 'lower .95']
        coef.upper <- model.summ$conf.int[var, 'upper .95']

        coef.summ.string <- .formatCoef(exp.coef, coef.lower, coef.upper)

        #Create a row
        row.output.df <- data.frame(var.name = var,
                                    coef = coef.summ.string,
                                    p.val = p.val)
        #Row names
        row.name <- ifelse(var.n == 1, model.name, '')
        rownames(row.output.df) <- row.name

        #Add the output row to the output DataFrame
        output.df <- rbind(output.df, row.output.df)

      }

    } else {

      #P value
      p.val <- PVal(model.summ$coefficients[independent.var, 'Pr(>|z|)'])

      #Get Exp(coef) and 95% CI from the model
      exp.coef <- model.summ$conf.int[independent.var, 'exp(coef)']
      coef.lower <- model.summ$conf.int[independent.var, 'lower .95']
      coef.upper <- model.summ$conf.int[independent.var, 'upper .95']

      if (styled) {
        coef.summ.string <- .formatCoef(exp.coef, coef.lower, coef.upper)

        #Create a row
        row.output.df <- data.frame(coef = coef.summ.string,
                                    p.val = p.val)

        #The rowname is either all variables together or from the model.names
        if (is.null(model.names) | length(model.names) != length(models)) {
        row.name <- toString(variables.vector)
        if (length(model.names) != length(models)) {
          warning('Length of model.names != length of models.')
        }
        } else {
          row.name <- model.names[model.n]
        }
        rownames(row.output.df) <- row.name

      } else {
        row.output.df <- data.frame(exp.coef = exp.coef,
                                    lower = coef.lower,
                                    upper = coef.upper,
                                    p.val = p.val)
      }

      #Add ouput row to the output DataFrame
      output.df <- rbind(output.df, row.output.df)

    }
  }

  return(output.df)

}

ForestPlotCox <- function(models,
                          label.vars  = NULL,
                          hrzl.lines = list('2' = grid::gpar(lty = 2)),
                          clip = NULL,
                          xticks = NULL,
                          graph.width = 'auto',
                          box.size = 0.2,
                          zero.loc = 1.0,
                          box.col = "#212427",
                          lines.col = "#212427",
                          zero.col = "red",
                          text.col = "#212427",
                          xlab = 'Hazard Ratios',
                          cex = 1.0,
                          cex.axis = 1.0,
                          cex.xlab = 1.0,
                          graph.pos = 2) {

  n.col.df <- data.frame()
  hr.df <- data.frame()
  labels.df <- data.frame()
  rows.n <- nrow(models)

  for (i in 1:rows.n) {
    n.col.df <- rbind(n.col.df, data.frame(model.name = paste0('Model ', i)))
  }

  for (i in 1:rows.n) {
    hr <- paste0(NumToString(models[i, "exp.coef"], 2), ' [', NumToString(models[i, "lower"], 2), ', ', NumToString(models[i, "upper"], 2), ']')
    hr.df <- rbind(hr.df, data.frame(hr))
  }

  if (is.null(label.vars) == FALSE) {
    for (label in label.vars) {
      labels.df <- rbind(labels.df, data.frame(label = paste0(label)))
    }
    models[, 'model.label'] <- labels.df
  } else {
    models[, 'model.label'] <- NA
  }

  models[, 'model.name'] <- n.col.df
  models[, 'hr'] <- hr.df
  names(models)[1] <- 'mean'
  models[, 'is.summary'] <- NA

  if (is.null(label.vars)) {
    header <- data.frame(model.name = 'Model',
                         model.label = NA,
                         p.val = 'P value',
                         hr = 'HR [95% CI]',
                         is.summary = TRUE,
                         mean = NA,
                         lower = NA,
                         upper = NA
    )
    label.text <- c('model.name', 'hr', 'p.val')
  } else {
    header <- data.frame(model.name = 'Model',
                         model.label = 'Variables in model',
                         hr = 'HR [95% CI]',
                         p.val = 'P value',
                         is.summary = TRUE,
                         mean = NA,
                         lower = NA,
                         upper = NA)
  }

  binded.df <- rbind(header, models)

  plot <- forestplot::forestplot(binded.df,
                                 labeltext = c(model.name, model.label, hr, p.val),
                                 xticks = xticks,
                                 clip = clip,
                                 graph.pos = graph.pos,
                                 zero = zero.loc,
                                 boxsize = box.size,
                                 hrzl_lines = hrzl.lines,
                                 is.summary = is.summary,
                                 xlab = xlab,
                                 col = forestplot::fpColors(box = box.col,
                                                lines = lines.col,
                                                zero = zero.col,
                                                text = text.col),
                                 graphwidth = graph.width,
                                 txt_gp = forestplot::fpTxtGp(label = grid::gpar(cex = cex),
                                                  ticks = grid::gpar(cex = cex.axis),
                                                  xlab = grid::gpar(cex = cex.xlab)))

  return(plot)

}

PlotHr <- function(pred,
                   df,
                   ind.var,
                   breaks = 15,
                   background.col = 'gray',
                   title.text = NULL,
                   round.x = 1,
                   round.y = 1,
                   axis.1.at = NULL,
                   axis.1.labels = NULL,
                   axis.2.at = NULL,
                   axis.2.labels = NULL,
                   axis.4.at = NULL,
                   axis.4.labels = NULL,
                   x.lab = 'x',
                   y.lab = 'y',
                   lim.x = NULL,
                   h.border = TRUE,
                   y.4.lab = 'Frequency (n)',
                   cex = 1.0,
                   cex.labels = 1.0,
                   mar.plot = c(5,5,2,7),
                   abline = TRUE,
                   abline.lty = 2,
                   box = TRUE,
                   abline.col = 'red',
                   polygon.col = 'black',
                   polygon.alpha = 0.3,
                   outline.polygon = FALSE,
                   h.wd = 1,
                   axes.wd = 1,
                   line.wd = 2) {

  #Define variables that will later help with plotting

  h <- hist(df[,colnames(df)==ind.var], breaks = breaks, plot = FALSE)

  ylim.h.top <- if (!is.null(axis.4.at)) {max(axis.4.at)} else {max(h$counts) + (max(h$counts) * 0.1)}

  ylim.bot <- min(pred$lower)
  ylim.top <- round(max(pred$upper))

  xlim.bot <- min(pred[, ind.var])
  xlim.top <- max(pred[, ind.var])

  if (is.null(lim.x)) {
    lim.x = c(xlim.bot, xlim.top)
  }

  if (is.null(axis.1.at)) {
    axis.1.at <- c(round(xlim.bot, round.x):round(xlim.top, round.x))
  }

  if (is.null(axis.1.labels)) {
    axis.1.labels <- c(round(xlim.bot, round.x):round(xlim.top, round.x))
  }

  if (is.null(axis.2.at)) {
    axis.2.at <- c(round(ylim.bot, round.y):round(ylim.top, round.y))
  }

  if (is.null(axis.2.labels)) {
    axis.2.labels <- c(round(ylim.bot, round.y):round(ylim.top, round.y))
  }

  dev.hold(); on.exit(dev.flush())

  #Create plot settings
  par(mar = mar.plot,
      xaxs = 'i',
      yaxs = 'i',
      cex.axis = cex)

  #Plot histogram behind
  plot.new()
  plot.window(xlim = lim.x, ylim = c(0,ylim.h.top))
  rect(h$breaks[-length(h$breaks)],  0, h$breaks[-1], h$counts,
       lwd = h.wd, border = h.border, angle = 45, lty = NULL, density = NULL)
  axis(1, at = axis.1.at, labels = axis.1.labels)
  axis(side = 4, if (!is.null(axis.4.at)) {at = axis.4.at})
  mtext(y.4.lab, side = 4, line = 3, las = 0, cex = cex.labels)
  if (is.null(title.text) == FALSE) {
    title(main = title.text)
  }

  #Plot the HR
  par(new = TRUE,
      cex.axis = cex)

  plot.window(xlim = lim.x,ylim = c(ylim.bot, ylim.top))
  axis(2, at = axis.2.at, labels = axis.2.labels)
  if (outline.polygon) {
    lines(pred[,colnames(pred)==ind.var], pred[,colnames(pred)=="upper"])
    lines(pred[,colnames(pred)==ind.var], pred[,colnames(pred)=="lower"])
  }
  polygon(x = c(pred[,colnames(pred)==ind.var], rev(pred[,colnames(pred)==ind.var])),
          y = c(pred[,colnames(pred)=="yhat"], rev(pred[,colnames(pred)=="upper"])),
          col = adjustcolor(polygon.col, polygon.alpha),
          border = NA)
  polygon(x = c(pred[,colnames(pred)==ind.var], rev(pred[,colnames(pred)==ind.var])),
          y = c(pred[,colnames(pred)=="yhat"], rev(pred[,colnames(pred)=="lower"])),
          col = adjustcolor(polygon.col, polygon.alpha),
          border = NA)
  lines(pred[,colnames(pred)==ind.var], pred[,colnames(pred)=="yhat"],lwd=line.wd)

  mtext(x.lab, side = 1, line = 3, las = 0, cex = cex.labels)
  mtext(y.lab, side = 2, line = 3, las = 0, cex = cex.labels)

  if (abline) {
    abline(h=1.0, col = abline.col, lty = abline.lty, lwd = line.wd)
  }

  if (box) {
    box()
  }

  invisible()

}

InteractCox <- function(independent.var,
                        adjust.for,
                        tests.interaction,
                        surv.obj,
                        df,
                        row.labels = NULL) {

  if (length(tests.interaction) != length(row.labels)) {
    warning('Length of row.labels does not much the length of tests.interaction.')
  }

  for (i in 1:length(adjust.for)) {
    if (i == 1) {
      adjust.for.formula <- adjust.for[i]
    } else if (i == length(adjust.for)) {
      adjust.for.formula <- paste0(adjust.for.formula,
                                   ' + ',
                                   adjust.for[i])
    } else {
      adjust.for.formula <- paste0(adjust.for.formula,
                                   ' + ',
                                   adjust.for[i],
                                   ' +')
    }
  }

  df.interactions <- data.frame(var = '',
                                numbers = 'N (N event)',
                                hr = 'HR [95% CI]',
                                p = 'P value for interaction',
                                is.summary = TRUE,
                                mean = NA,
                                lower = NA,
                                upper = NA)
  iteration <- 0

  for (test.interaction in tests.interaction) {

    iteration <- iteration + 1

    fit.1.formula <- as.formula(paste0(deparse(substitute(surv.obj)),
                                       ' ~ ',
                                       independent.var,
                                       ' + ',
                                       adjust.for.formula,
                                       ' + ',
                                       test.interaction[1]
    ))
    fit.1 <- survival::coxph(fit.1.formula, data = df)

    fit.2.formula <- as.formula(paste0(deparse(substitute(surv.obj)),
                                       ' ~ ',
                                       independent.var,
                                       ' + ',
                                       adjust.for.formula,
                                       ' + ',
                                       test.interaction[1],
                                       ' + ',
                                       independent.var,
                                       ':',
                                       test.interaction[1]
    ))
    fit.2 <- survival::coxph(fit.2.formula, data = df)
    p.interaction <- round(anova(fit.1, fit.2)$'P(>|Chi|)'[2], 3)

    fit.l.formula <- as.formula(paste0(deparse(substitute(surv.obj)),
                                       ' ~ ',
                                       independent.var,
                                       ' + ',
                                       adjust.for.formula,
                                       ' '))
    fit.l <- survival::coxph(fit.l.formula,
                             subset = (df[,test.interaction[1]]  == test.interaction[2]),
                             data = df)
    fit.l.hr <- c(round(summary(fit.l)$conf.int[independent.var, 'exp(coef)'], 2),
                  round(summary(fit.l)$conf.int[independent.var, 'lower .95'], 2),
                  round(summary(fit.l)$conf.int[independent.var, 'upper .95'], 2)
    )

    fit.h.formula <- as.formula(paste0(deparse(substitute(surv.obj)),
                                       ' ~ ',
                                       independent.var,
                                       ' + ',
                                       adjust.for.formula,
                                       ' '))
    fit.h <- survival::coxph(fit.l.formula,
                             subset = (df[,test.interaction[1]]  == test.interaction[3]),
                             data = df)
    fit.h.hr <- c(round(summary(fit.h)$conf.int[independent.var, 'exp(coef)'], 2),
                  round(summary(fit.h)$conf.int[independent.var, 'lower .95'], 2),
                  round(summary(fit.h)$conf.int[independent.var, 'upper .95'], 2)
    )

    #Check what to use for the header of the rows
    name.rows <- if (is.na(row.labels[iteration]) | is.null(row.labels[iteration])) {
      test.interaction[1]
    } else {
      row.labels[iteration]
    }

    row <- data.frame(var = c(name.rows,
                              test.interaction[2],
                              test.interaction[3]),
                      numbers = c(NA,
                                  paste0(summary(fit.l)$n,
                                         ' (',
                                         summary(fit.l)$nevent,
                                         ')'),
                                  paste0(summary(fit.h)$n,
                                         ' (',
                                         summary(fit.h)$nevent,
                                         ')')),
                      hr = c(NA,
                             paste0(fit.l.hr[1],
                                    ' [',
                                    fit.l.hr[2],
                                    ', ',
                                    fit.l.hr[3],
                                    ']'),
                             paste0(fit.h.hr[1],
                                    ' [',
                                    fit.h.hr[2],
                                    ', ',
                                    fit.h.hr[3],
                                    ']')
                      ),
                      p = c(p.interaction, NA, NA),
                      is.summary = c(TRUE, NA, NA),
                      mean = c(NA, fit.l.hr[1],fit.h.hr[1]),
                      lower = c(NA, fit.l.hr[2],fit.h.hr[2]),
                      upper = c(NA, fit.l.hr[3],fit.h.hr[3]))

    df.interactions <- rbind(df.interactions, row)
  }

  return(df.interactions)
}



