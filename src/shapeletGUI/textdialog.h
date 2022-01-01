#ifndef TEXTDIALOG_H
#define TEXTDIALOG_H

#include <QDialog>

namespace Ui {
class TextDialog;
}

class TextDialog : public QDialog
{
    Q_OBJECT

public:
    explicit TextDialog(QWidget *parent = nullptr,  const QString& msg="");
    ~TextDialog();

private:
    Ui::TextDialog *ui;
};

#endif // TEXTDIALOG_H
