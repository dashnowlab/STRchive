import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";

export default function ContactLink() {
  return (
    <Dialog
      title="Contact Us"
      trigger={
        <Button id="footer-contact" className="p-0! underline">
          Contact
        </Button>
      }
    >
      <ContactForm />
    </Dialog>
  );
}
